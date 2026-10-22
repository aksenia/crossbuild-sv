#!/usr/bin/env python3
"""
Collapse sorted, normalized VCF records with identical CHROM/POS/END/ID,
joining their LIFTOVER_TOOL values into a comma-separated list.

The ID was stamped as CHROM_POS_SVTYPE_END in hg19 coordinates before
liftover, so it identifies the source variant; CHROM/POS/END is where a tool
put it on hg38. Inputs are normalized first (normalize_liftover), so equal
coordinates mean the tools agree and any difference is a real disagreement.
REF/ALT are not compared (SV alleles can be kb long); the first tool in
TOOL_ORDER supplies the record.

Also:
- Carries tool-specific INFO fields (e.g. transanno MULTIMAP/REF_CHANGED).
- Optionally annotates records with Picard rejection metadata when
  --picard-rejected is supplied.
- With --refchanged-ids, every record of a listed source ID (the reference
  allele differs between builds, flagged by at least one tool) goes to
  --refchanged-out instead, tagged REFCHANGED_BY=<flagging tools>. Routing is
  per ID across all tools, so a tool that silently mis-handles such a locus
  does not leave its record in the main set.

Reads from stdin (sorted VCF), writes to stdout.
"""
import argparse
import sys

TOOL_ORDER = ('crossmap', 'bcftools', 'picard', 'transanno')

_TOOL_FIELDS = {
    'transanno': ('MULTIMAP', 'REF_CHANGED'),
}

_PICARD_HEADERS = [
    '##INFO=<ID=PICARD_REJECT_REASON,Number=1,Type=String,Description="Picard LiftoverVcf rejection reason">',
    '##INFO=<ID=PICARD_ATTEMPTED_LOCUS,Number=1,Type=String,Description="hg38 locus Picard attempted before failing (MismatchedRefAllele only)">',
    '##INFO=<ID=PICARD_ATTEMPTED_ALLELES,Number=1,Type=String,Description="hg38 alleles Picard attempted before failing (MismatchedRefAllele only)">',
]

_REFCHANGED_HEADER = '##INFO=<ID=REFCHANGED_BY,Number=.,Type=String,Description="Tools that flagged this source variant as REF-changed between builds (multiallelic, SWAP or REF_CHANGED)">'


def parse_info(info_str):
    d = {}
    for field in info_str.split(';'):
        if '=' in field:
            k, v = field.split('=', 1)
            d[k] = v
        elif field:
            d[field] = True
    return d


def format_info(d):
    return ';'.join(k if v is True else f'{k}={v}' for k, v in d.items())


def tools_of(info):
    return [t for t in str(info.get('LIFTOVER_TOOL', '')).split(',') if t]


def tool_rank(tool):
    return TOOL_ORDER.index(tool) if tool in TOOL_ORDER else len(TOOL_ORDER)


def load_picard_rejected(path):
    """Parse Picard rejected VCF → {variant_id: {reason, locus, alleles}}."""
    rejected = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 8:
                continue
            info = parse_info(parts[7])
            rejected[parts[2]] = {
                'reason':  parts[6],
                'locus':   info.get('AttemptedLocus', ''),
                'alleles': info.get('AttemptedAlleles', ''),
            }
    return rejected


def _apply_picard(variant_id, info, picard_rejected):
    rec = picard_rejected.get(variant_id)
    if not rec:
        return
    info['PICARD_REJECT_REASON'] = rec['reason']
    if rec['locus']:
        info['PICARD_ATTEMPTED_LOCUS'] = rec['locus']
    if rec['alleles']:
        info['PICARD_ATTEMPTED_ALLELES'] = rec['alleles']


def load_refchanged(path):
    """'ID<TAB>tool' lines -> {ID: {tools}}."""
    flagged = {}
    with open(path) as fh:
        for line in fh:
            vid, tool = line.split()
            flagged.setdefault(vid, set()).add(tool)
    return flagged


def flush_group(records, picard_rejected, refchanged, rc_out):
    """Write one merged record for records sharing CHROM/POS/END/ID."""
    infos = [parse_info(p[7]) for p in records]
    order = sorted(range(len(records)), key=lambda i: min(map(tool_rank, tools_of(infos[i])), default=99))
    merged, info = list(records[order[0]]), infos[order[0]]

    tools = {t for i in infos for t in tools_of(i)}
    for rec_info in infos:
        for tool_name, fields in _TOOL_FIELDS.items():
            if tool_name in tools_of(rec_info):
                info.update({f: rec_info[f] for f in fields if f in rec_info})

    info['LIFTOVER_TOOL'] = ','.join(sorted(tools, key=tool_rank))
    _apply_picard(merged[2], info, picard_rejected)
    out = sys.stdout
    if merged[2] in refchanged:
        info['REFCHANGED_BY'] = ','.join(sorted(refchanged[merged[2]], key=tool_rank))
        out = rc_out
    merged[7] = format_info(info)
    out.write('\t'.join(merged) + '\n')


def flush_position(records, *args):
    """
    Group all records at one CHROM/POS by END/ID. Sorting is by CHROM/POS
    only, so records of different IDs can interleave at the same position;
    grouping consecutive lines alone would split one ID into several records.
    """
    groups = {}
    for parts in records:
        end = parse_info(parts[7])['END']   # always set by finalize_sv.py
        groups.setdefault((int(end), parts[2]), []).append(parts)
    for key in sorted(groups):
        flush_group(groups[key], *args)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--picard-rejected', metavar='VCF',
                        help='Picard rejected VCF; annotates records with rejection reason')
    parser.add_argument('--refchanged-ids', metavar='TSV',
                        help="'ID<TAB>tool' lines for source variants flagged REF-changed")
    parser.add_argument('--refchanged-out', metavar='VCF',
                        help='plain VCF receiving records of --refchanged-ids')
    args = parser.parse_args()
    if bool(args.refchanged_ids) != bool(args.refchanged_out):
        parser.error('--refchanged-ids and --refchanged-out go together')

    picard_rejected = load_picard_rejected(args.picard_rejected) if args.picard_rejected else {}
    picard_headers_injected = not picard_rejected
    refchanged = load_refchanged(args.refchanged_ids) if args.refchanged_ids else {}
    rc_out = open(args.refchanged_out, 'w') if args.refchanged_out else None
    state = (picard_rejected, refchanged, rc_out)

    prev_pos = None
    bucket = []

    for raw_line in sys.stdin:
        line = raw_line.rstrip('\n')
        if line.startswith('#'):
            if not picard_headers_injected and line.startswith('#CHROM'):
                for h in _PICARD_HEADERS:
                    sys.stdout.write(h + '\n')
                    if rc_out:
                        rc_out.write(h + '\n')
                picard_headers_injected = True
            sys.stdout.write(line + '\n')
            if rc_out:
                if line.startswith('#CHROM'):
                    rc_out.write(_REFCHANGED_HEADER + '\n')
                rc_out.write(line + '\n')
            continue

        parts = line.split('\t')
        pos = (parts[0], parts[1])
        if pos != prev_pos and bucket:
            flush_position(bucket, *state)
            bucket = []
        bucket.append(parts)
        prev_pos = pos

    if bucket:
        flush_position(bucket, *state)
    if rc_out:
        rc_out.close()


if __name__ == '__main__':
    main()
