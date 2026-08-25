#!/usr/bin/env python3
"""Normalize only the target END field in lifted structural-variant VCFs.

The liftover tool's CHROM, POS, ID, REF and ALT are never rewritten.

Sequence-resolved records:
    END = POS + len(REF) - 1

Symbolic INS records:
    END = POS

Symbolic DEL/DUP/DUP:TANDEM/INV records:
    Map source POS and END through one chain record.  The mapped interval must
    have the same target contig and POS already reported by the liftover tool;
    only then is END replaced by the mapped target endpoint.

Records that cannot satisfy these rules are written to a separate rejection
VCF instead of silently retaining a source-build END.
"""

import argparse
import bisect
import gzip
import re
import sys
from dataclasses import dataclass

import pysam


INTERVAL_TYPES = {"DEL", "DUP", "DUP:TANDEM", "INV"}


def open_text(path, mode="rt"):
    path = str(path)
    if path.endswith((".gz", ".bgz")):
        return gzip.open(path, mode)
    return open(path, mode)


def parse_info(text):
    info = {}
    if text in ("", "."):
        return info
    for item in text.split(";"):
        if not item:
            continue
        if "=" in item:
            key, value = item.split("=", 1)
            info[key] = value
        else:
            info[item] = True
    return info


def format_info(info):
    if not info:
        return "."
    fields = []
    for key, value in info.items():
        if value is True:
            fields.append(key)
        elif value is not None:
            fields.append(f"{key}={value}")
    return ";".join(fields) if fields else "."


def first_int(value):
    if value in (None, "", "."):
        return None
    try:
        return int(str(value).split(",", 1)[0])
    except ValueError:
        return None


def bare_contig(chrom):
    if chrom == "chrM":
        return "MT"
    return chrom[3:] if chrom.startswith("chr") else chrom


def equivalent_contig(first, second):
    return first == second or bare_contig(first) == bare_contig(second)


def is_symbolic_alt(alt):
    return alt.startswith("<") and alt.endswith(">")


def representation(alts):
    values = alts.split(",") if alts else []
    if not values:
        return "unsupported"
    if any(alt in (".", "*") or "[" in alt or "]" in alt for alt in values):
        return "unsupported"
    symbolic = [is_symbolic_alt(alt) for alt in values]
    if all(symbolic):
        return "symbolic"
    if not any(symbolic):
        return "sequence"
    return "unsupported"


@dataclass
class SourceVariant:
    chrom: str
    pos: int
    end: int | None
    svtype: str
    representation: str


def load_source_variants(path):
    variants = {}
    with open_text(path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8:
                continue

            chrom, pos_text, variant_id, _ref, alt = parts[:5]
            if variant_id in ("", "."):
                raise ValueError("Source VCF contains a record without a stamped ID")
            if variant_id in variants:
                raise ValueError(f"Duplicate source ID after preprocessing: {variant_id}")

            info = parse_info(parts[7])
            variants[variant_id] = SourceVariant(
                chrom=chrom,
                pos=int(pos_text),
                end=first_int(info.get("END")),
                svtype=str(info.get("SVTYPE", "")),
                representation=representation(alt),
            )
    return variants


@dataclass
class ChainBlock:
    source_start: int
    source_end: int
    target_start: int


@dataclass
class ChainRecord:
    score: int
    chain_id: str
    source_name: str
    source_strand: str
    target_name: str
    target_size: int
    target_strand: str
    blocks: list

    def __post_init__(self):
        self.refresh()

    def refresh(self):
        self._starts = [block.source_start for block in self.blocks]

    def map_point(self, position):
        """Map one 1-based source base; return a 1-based target base."""
        if self.source_strand != "+":
            return None

        source_zero = position - 1
        index = bisect.bisect_right(self._starts, source_zero) - 1
        if index < 0:
            return None

        block = self.blocks[index]
        if not (block.source_start <= source_zero < block.source_end):
            return None

        offset = source_zero - block.source_start
        target_oriented = block.target_start + offset
        if self.target_strand == "+":
            target_zero = target_oriented
        else:
            target_zero = self.target_size - target_oriented - 1
        return target_zero + 1


class ChainIndex:
    def __init__(self, path):
        self.by_source = {}
        self._load(path)

    def _store(self, record):
        record.refresh()
        self.by_source.setdefault(record.source_name, []).append(record)

    def _load(self, path):
        current = None
        source_position = target_position = None

        with open_text(path) as handle:
            for raw in handle:
                line = raw.strip()
                if not line:
                    if current is not None:
                        self._store(current)
                        current = None
                    continue

                fields = line.split()
                if fields[0] == "chain":
                    if current is not None:
                        self._store(current)
                    if len(fields) != 13:
                        raise ValueError(f"Invalid chain header: {line}")

                    (
                        _, score, source_name, _source_size, source_strand,
                        source_start, _source_end, target_name, target_size,
                        target_strand, target_start, _target_end, chain_id,
                    ) = fields
                    current = ChainRecord(
                        score=int(score),
                        chain_id=chain_id,
                        source_name=source_name,
                        source_strand=source_strand,
                        target_name=target_name,
                        target_size=int(target_size),
                        target_strand=target_strand,
                        blocks=[],
                    )
                    source_position = int(source_start)
                    target_position = int(target_start)
                    continue

                if current is None:
                    raise ValueError("Chain block encountered before chain header")

                size = int(fields[0])
                current.blocks.append(
                    ChainBlock(
                        source_start=source_position,
                        source_end=source_position + size,
                        target_start=target_position,
                    )
                )
                source_position += size
                target_position += size

                if len(fields) == 3:
                    source_position += int(fields[1])
                    target_position += int(fields[2])
                elif len(fields) != 1:
                    raise ValueError(f"Invalid chain block line: {line}")

        if current is not None:
            self._store(current)

    def source_records(self, chrom):
        aliases = [chrom]
        if chrom.startswith("chr"):
            aliases.append(bare_contig(chrom))
        else:
            aliases.append("chrM" if chrom == "MT" else f"chr{chrom}")

        records = []
        seen = set()
        for alias in aliases:
            for record in self.by_source.get(alias, []):
                key = (record.chain_id, record.source_name, record.target_name)
                if key not in seen:
                    seen.add(key)
                    records.append(record)
        return records

    def target_end_for_tool_position(
        self, source_chrom, source_pos, source_end, target_chrom, target_pos
    ):
        """Return END only when a same-chain pair agrees with tool CHROM/POS."""
        matches = []
        for chain in self.source_records(source_chrom):
            mapped_pos = chain.map_point(source_pos)
            mapped_end = chain.map_point(source_end)
            if mapped_pos is None or mapped_end is None:
                continue
            if not equivalent_contig(chain.target_name, target_chrom):
                continue

            interval_pos = min(mapped_pos, mapped_end)
            interval_end = max(mapped_pos, mapped_end)
            if interval_pos != target_pos:
                continue
            matches.append((-chain.score, chain.chain_id, interval_end))

        if not matches:
            return None
        matches.sort()
        return matches[0][2]


class FastaLookup:
    def __init__(self, path):
        self.fasta = pysam.FastaFile(path)
        self.names = set(self.fasta.references)

    def resolve_contig(self, chrom):
        if chrom in self.names:
            return chrom
        bare = bare_contig(chrom)
        candidates = [bare, "chrM" if bare == "MT" else f"chr{bare}"]
        for candidate in candidates:
            if candidate in self.names:
                return candidate
        return None

    def validate_ref(self, chrom, position, ref):
        resolved = self.resolve_contig(chrom)
        if resolved is None:
            raise ValueError("TARGET_CONTIG_NOT_IN_FASTA")
        observed = self.fasta.fetch(
            resolved, position - 1, position - 1 + len(ref)
        ).upper()
        if len(observed) != len(ref):
            raise ValueError("TARGET_REF_OUT_OF_RANGE")
        if observed != ref.upper():
            raise ValueError("TARGET_REF_MISMATCH")


def inject_headers(header_lines):
    additions = [
        '##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant on the target reference">',
        '##INFO=<ID=END_LIFTOVER_REASON,Number=1,Type=String,Description="Reason a lifted record failed END normalization">',
        '##FILTER=<ID=END_LIFTOVER_FAILED,Description="Target SV END could not be established without changing the tool result">',
    ]
    existing = "\n".join(header_lines)
    insert_at = next(
        (index for index, line in enumerate(header_lines) if line.startswith("#CHROM")),
        len(header_lines),
    )
    for addition in additions:
        token = addition.split("<ID=", 1)[1].split(",", 1)[0]
        if f"<ID={token}," not in existing:
            header_lines.insert(insert_at, addition)
            insert_at += 1
            existing += "\n" + addition
    return header_lines


def set_filter(current, new_filter):
    if current in ("", ".", "PASS"):
        return new_filter
    values = current.split(";")
    if new_filter not in values:
        values.append(new_filter)
    return ";".join(values)


def safe_reason(text):
    cleaned = re.sub(r"[^A-Za-z0-9_.:-]+", "_", text)
    return cleaned[:200] or "UNKNOWN"


def normalize_record(parts, source, chains, fasta):
    chrom = parts[0]
    position = int(parts[1])
    ref = parts[3]
    info = parse_info(parts[7])

    # Validation is non-mutating: the tool's target REF is either retained as
    # reported or the record is rejected.
    fasta.validate_ref(chrom, position, ref)

    if source.representation == "sequence":
        target_end = position + len(ref) - 1
    elif source.representation == "symbolic" and source.svtype == "INS":
        target_end = position
    elif source.representation == "symbolic" and source.svtype in INTERVAL_TYPES:
        if source.end is None:
            raise ValueError("SOURCE_END_UNAVAILABLE")
        target_end = chains.target_end_for_tool_position(
            source.chrom,
            source.pos,
            source.end,
            chrom,
            position,
        )
        if target_end is None:
            raise ValueError("CHAIN_PAIR_DOES_NOT_AGREE_WITH_TOOL_POSITION")
    else:
        raise ValueError("UNSUPPORTED_SV_REPRESENTATION")

    if target_end < position:
        raise ValueError("END_LT_POS_AFTER_NORMALIZATION")

    info["END"] = str(target_end)
    parts[7] = format_info(info)
    return parts


def read_lifted(path):
    headers = []
    records = []
    with open_text(path) as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith("#"):
                headers.append(line)
            elif line:
                records.append(line)
    return inject_headers(headers), records


def rejected_record(parts, reason):
    failed = parts[:]
    info = parse_info(failed[7])
    info["END_LIFTOVER_REASON"] = safe_reason(reason)
    failed[6] = set_filter(failed[6], "END_LIFTOVER_FAILED")
    failed[7] = format_info(info)
    return failed


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", required=True, help="ID-stamped source hg19 VCF")
    parser.add_argument("--lifted", required=True, help="Raw VCF from one liftover tool")
    parser.add_argument("--chain", required=True, help="Chain used by that tool branch")
    parser.add_argument("--target-fasta", required=True, help="Target hg38 FASTA")
    parser.add_argument("--tool", required=True)
    parser.add_argument("--output", required=True, help="END-normalized successful VCF")
    parser.add_argument("--reject", required=True, help="END-normalization rejection VCF")
    args = parser.parse_args()

    source_variants = load_source_variants(args.source)
    chains = ChainIndex(args.chain)
    fasta = FastaLookup(args.target_fasta)
    headers, records = read_lifted(args.lifted)

    accepted = 0
    rejected = 0
    with open(args.output, "w") as output, open(args.reject, "w") as reject:
        for header in headers:
            output.write(header + "\n")
            reject.write(header + "\n")

        for raw in records:
            parts = raw.split("\t")
            if len(parts) < 8:
                continue

            source = source_variants.get(parts[2])
            if source is None:
                reject.write(
                    "\t".join(rejected_record(parts, "SOURCE_ID_NOT_FOUND")) + "\n"
                )
                rejected += 1
                continue

            try:
                normalized = normalize_record(parts[:], source, chains, fasta)
            except Exception as error:
                reject.write("\t".join(rejected_record(parts, str(error))) + "\n")
                rejected += 1
                continue

            output.write("\t".join(normalized) + "\n")
            accepted += 1

    print(
        f"[{args.tool}] END normalization: {accepted} accepted, "
        f"{rejected} rejected",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
