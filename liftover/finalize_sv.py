#!/usr/bin/env python3
"""
Finalize one normalized SV VCF (one liftover tool's output, or native hg38).

Input must already be through `bcftools norm -f <hg38>` (left-aligned,
shared padding trimmed), so every tool's alleles are in the same canonical
representation. This script then:

- validates REF against the hg38 FASTA (mismatches are rejected, not repaired);
- sets INFO/END from the record itself:
    sequence-resolved   END = POS + len(REF) - 1
    symbolic INS        END = POS
    symbolic DEL/DUP/INV  source POS and END mapped through one chain record;
                        accepted only if the mapped start equals the tool's POS
                        (native mode: the record's own END is kept)
- tags INFO/LIFTOVER_TOOL=<tool> (tool mode only).

Reads a VCF from stdin, writes plain VCFs.
"""
import argparse
import bisect
import gzip
import re
import sys
from dataclasses import dataclass

import pysam

INTERVAL_TYPES = {'DEL', 'DUP', 'INV'}

_HEADERS = {
    'END': '##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant on hg38">',
    'LIFTOVER_TOOL': '##INFO=<ID=LIFTOVER_TOOL,Number=.,Type=String,Description="Liftover tool(s) that produced this record">',
    'QC_REJECT_REASON': '##INFO=<ID=QC_REJECT_REASON,Number=1,Type=String,Description="Why the record failed post-liftover QC">',
    'QC_FAILED': '##FILTER=<ID=QC_FAILED,Description="Failed post-liftover QC (REF or END could not be established)">',
}


def open_text(path):
    return gzip.open(path, 'rt') if str(path).endswith('.gz') else open(path)


def parse_info(info_str):
    d = {}
    if info_str in ('', '.'):
        return d
    for field in info_str.split(';'):
        if '=' in field:
            k, v = field.split('=', 1)
            d[k] = v
        elif field:
            d[field] = True
    return d


def format_info(d):
    if not d:
        return '.'
    return ';'.join(k if v is True else f'{k}={v}' for k, v in d.items())


def first_int(value):
    try:
        return int(str(value).split(',', 1)[0])
    except (TypeError, ValueError):
        return None


def bare_contig(chrom):
    if chrom == 'chrM':
        return 'MT'
    return chrom[3:] if chrom.startswith('chr') else chrom


def is_symbolic(alt):
    return alt.startswith('<') and alt.endswith('>')


@dataclass
class SourceVariant:
    chrom: str
    pos: int
    end: int


def load_source(path):
    """Stamped hg19 source VCF -> {ID: SourceVariant}."""
    variants = {}
    with open_text(path) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            info = parse_info(parts[7])
            end = first_int(info.get('END'))
            variants[parts[2]] = SourceVariant(parts[0], int(parts[1]),
                                               end if end is not None else int(parts[1]))
    return variants


@dataclass
class ChainBlock:
    source_start: int
    source_end: int
    target_start: int


class Chain:
    def __init__(self, score, target_name, target_size, target_strand):
        self.score = score
        self.target_name = target_name
        self.target_size = target_size
        self.target_strand = target_strand
        self.blocks = []
        self.starts = []

    def map_point(self, pos):
        """Map a 1-based source base to a 1-based target base, or None."""
        zero = pos - 1
        i = bisect.bisect_right(self.starts, zero) - 1
        if i < 0:
            return None
        block = self.blocks[i]
        if not block.source_start <= zero < block.source_end:
            return None
        target = block.target_start + zero - block.source_start
        if self.target_strand == '-':
            target = self.target_size - target - 1
        return target + 1


class ChainIndex:
    """UCSC-format chain file indexed by source contig (+ strand sources only)."""

    def __init__(self, path):
        self.by_source = {}
        chain = None
        with open_text(path) as fh:
            for line in fh:
                fields = line.split()
                if not fields:
                    continue
                if fields[0] == 'chain':
                    chain = None
                    if fields[4] != '+':
                        continue
                    chain = Chain(int(fields[1]), fields[7], int(fields[8]), fields[9])
                    s, t = int(fields[5]), int(fields[10])
                    self.by_source.setdefault(bare_contig(fields[2]), []).append(chain)
                    continue
                if chain is None:
                    continue
                size = int(fields[0])
                chain.blocks.append(ChainBlock(s, s + size, t))
                chain.starts.append(s)
                s += size
                t += size
                if len(fields) == 3:
                    s += int(fields[1])
                    t += int(fields[2])

    def target_end(self, src, target_chrom, target_pos):
        """END for a symbolic interval, if one chain maps it onto the tool's POS."""
        best = None
        for chain in self.by_source.get(bare_contig(src.chrom), []):
            if bare_contig(chain.target_name) != bare_contig(target_chrom):
                continue
            a, b = chain.map_point(src.pos), chain.map_point(src.end)
            if a is None or b is None or min(a, b) != target_pos:
                continue
            if best is None or chain.score > best[0]:
                best = (chain.score, max(a, b))
        return best[1] if best else None


def target_end(parts, info, source, chains):
    pos, ref, alts = int(parts[1]), parts[3], parts[4].split(',')
    svtype = str(info.get('SVTYPE', ''))
    if any(a in ('.', '*') or '[' in a or ']' in a for a in alts):
        raise ValueError('UNSUPPORTED_ALT')
    symbolic = [is_symbolic(a) for a in alts]
    if not any(symbolic):
        return pos + len(ref) - 1
    if not all(symbolic):
        raise ValueError('MIXED_SYMBOLIC_ALT')
    if svtype == 'INS':
        return pos
    if svtype not in INTERVAL_TYPES:
        raise ValueError('UNSUPPORTED_SVTYPE')
    if chains is None:
        end = first_int(info.get('END'))
        if end is None:
            raise ValueError('END_MISSING')
        return end
    src = source.get(parts[2])
    if src is None:
        raise ValueError('SOURCE_ID_NOT_FOUND')
    end = chains.target_end(src, parts[0], pos)
    if end is None:
        raise ValueError('CHAIN_DOES_NOT_MATCH_TOOL_POS')
    return end


def check_ref(fasta, chrom, pos, ref):
    if chrom not in fasta.references:
        raise ValueError('CONTIG_NOT_IN_HG38')
    observed = fasta.fetch(chrom, pos - 1, pos - 1 + len(ref)).upper()
    if observed != ref.upper():
        raise ValueError('REF_MISMATCH')


def reject(parts, info, reason):
    info['QC_REJECT_REASON'] = re.sub(r'[^A-Za-z0-9_.:-]+', '_', reason)[:200]
    parts[6] = 'QC_FAILED'
    parts[7] = format_info(info)
    return parts


def with_headers(header, tool):
    wanted = ['END', 'QC_REJECT_REASON', 'QC_FAILED'] + (['LIFTOVER_TOOL'] if tool else [])
    present = '\n'.join(header)
    extra = [_HEADERS[k] for k in wanted if f'<ID={k},' not in present]
    if tool and '<ID=LIFTOVER_TOOL,' in present:
        header = [h for h in header if '<ID=LIFTOVER_TOOL,' not in h]
        extra.append(_HEADERS['LIFTOVER_TOOL'])
    return header[:-1] + extra + header[-1:]


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('--tool', help='Liftover tool name; omit for native hg38 calls')
    parser.add_argument('--source', help='Stamped hg19 source VCF (tool mode)')
    parser.add_argument('--chain', help='Chain used by this tool (tool mode)')
    parser.add_argument('--fasta', required=True, help='hg38 FASTA')
    parser.add_argument('--output', required=True)
    parser.add_argument('--reject', required=True)
    args = parser.parse_args()
    if args.tool and not (args.source and args.chain):
        parser.error('--tool requires --source and --chain')

    source = load_source(args.source) if args.tool else {}
    chains = ChainIndex(args.chain) if args.tool else None
    fasta = pysam.FastaFile(args.fasta)

    counts = {'kept': 0, 'rejected': 0}
    header = []
    with open(args.output, 'w') as out, open(args.reject, 'w') as rej:
        for line in sys.stdin:
            if line.startswith('#'):
                header.append(line.rstrip('\n'))
                if line.startswith('#CHROM'):
                    for fh in (out, rej):
                        fh.write('\n'.join(with_headers(header, args.tool)) + '\n')
                continue
            parts = line.rstrip('\n').split('\t')
            info = parse_info(parts[7])
            try:
                check_ref(fasta, parts[0], int(parts[1]), parts[3])
                end = target_end(parts, info, source, chains)
                if end < int(parts[1]):
                    raise ValueError('END_BEFORE_POS')
            except ValueError as err:
                rej.write('\t'.join(reject(parts, info, str(err))) + '\n')
                counts['rejected'] += 1
                continue
            info['END'] = str(end)
            if args.tool:
                info['LIFTOVER_TOOL'] = args.tool
            parts[7] = format_info(info)
            out.write('\t'.join(parts) + '\n')
            counts['kept'] += 1

    print(f"[{args.tool or 'native'}] finalize: " +
          ', '.join(f'{k}={v}' for k, v in counts.items()), file=sys.stderr)


if __name__ == '__main__':
    main()
