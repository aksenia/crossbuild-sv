#!/usr/bin/env python3
"""
Merge coordinate-sorted liftover VCF records.

Records are collapsed only when they agree on:
    CHROM / POS / END / source-derived ID

This preserves SV interval disagreements while keeping the original pipeline's
coordinate-focused merge.  REF and ALT are not agreement criteria; the full
per-tool records remain available in the per-tool VCFs.

Input must be coordinate sorted.  The script buffers one CHROM/POS at a time,
then groups records within that position by the full merge key.
"""

import argparse
import gzip
import sys
from collections import OrderedDict


_TOOL_FIELDS = {
    "transanno": ("MULTIMAP", "REF_CHANGED"),
}

_TOOL_ORDER = ("crossmap", "bcftools", "picard", "transanno")

_PICARD_HEADERS = [
    '##INFO=<ID=PICARD_REJECT_REASON,Number=1,Type=String,Description="Picard LiftoverVcf rejection reason">',
    '##INFO=<ID=PICARD_ATTEMPTED_LOCUS,Number=1,Type=String,Description="hg38 locus Picard attempted before failing (MismatchedRefAllele only)">',
    '##INFO=<ID=PICARD_ATTEMPTED_ALLELES,Number=1,Type=String,Description="hg38 alleles Picard attempted before failing (MismatchedRefAllele only)">',
]


def open_text(path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


def parse_info(info_str):
    d = OrderedDict()

    if info_str in ("", "."):
        return d

    for field in info_str.split(";"):
        if not field:
            continue
        if "=" in field:
            k, v = field.split("=", 1)
            d[k] = v
        else:
            d[field] = True

    return d


def format_info(d):
    if not d:
        return "."

    parts = []
    for k, v in d.items():
        if v is True:
            parts.append(k)
        else:
            parts.append(f"{k}={v}")

    return ";".join(parts)


def record_end(parts):
    info = parse_info(parts[7])
    value = info.get("END")

    if value not in (None, "", "."):
        try:
            return int(str(value).split(",", 1)[0])
        except ValueError:
            pass

    raise ValueError(
        f"Normalized record lacks an integer END: {parts[0]}:{parts[1]}:{parts[2]}"
    )


def merge_key(parts):
    return (
        parts[0],
        int(parts[1]),
        record_end(parts),
        parts[2],
    )


def load_picard_rejected(path):
    """Parse Picard rejected VCF -> {variant_id: {reason, locus, alleles}}."""
    rejected = {}

    with open_text(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue

            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8:
                continue

            vid = parts[2]
            reason = parts[6]
            info = parse_info(parts[7])

            rejected[vid] = {
                "reason": reason,
                "locus": info.get("AttemptedLocus", ""),
                "alleles": info.get("AttemptedAlleles", ""),
            }

    return rejected


def _apply_picard(variant_id, info, picard_rejected):
    if not picard_rejected or variant_id not in picard_rejected:
        return

    rec = picard_rejected[variant_id]
    info["PICARD_REJECT_REASON"] = rec["reason"]

    if rec["locus"]:
        info["PICARD_ATTEMPTED_LOCUS"] = rec["locus"]

    if rec["alleles"]:
        info["PICARD_ATTEMPTED_ALLELES"] = rec["alleles"]


def ordered_tools(tool_names):
    seen = set(tool_names)
    result = [tool for tool in _TOOL_ORDER if tool in seen]
    result.extend(sorted(seen.difference(_TOOL_ORDER)))
    return result


def representative_rank(parts):
    """Choose a deterministic payload when tools agree on coordinates."""
    info = parse_info(parts[7])
    tools = ordered_tools(
        tool.strip()
        for tool in str(info.get("LIFTOVER_TOOL", "")).split(",")
        if tool.strip()
    )
    rank = (
        _TOOL_ORDER.index(tools[0])
        if tools and tools[0] in _TOOL_ORDER
        else len(_TOOL_ORDER)
    )
    return (rank, parts[3], parts[4], "\t".join(parts))


def flush_group(records, picard_rejected):
    if not records:
        return

    merged = list(min(records, key=representative_rank))
    merged_info = parse_info(merged[7])

    tools = []

    for parts in records:
        rec_info = parse_info(parts[7])

        for tool in str(rec_info.get("LIFTOVER_TOOL", "")).split(","):
            tool = tool.strip()
            if tool:
                tools.append(tool)

        tool_string = str(rec_info.get("LIFTOVER_TOOL", ""))

        for tool_name, fields in _TOOL_FIELDS.items():
            if tool_name in tool_string:
                for field in fields:
                    if field in rec_info:
                        merged_info[field] = rec_info[field]

    tool_union = ordered_tools(tools)
    if tool_union:
        merged_info["LIFTOVER_TOOL"] = ",".join(tool_union)

    _apply_picard(merged[2], merged_info, picard_rejected)

    merged[7] = format_info(merged_info)
    sys.stdout.write("\t".join(merged) + "\n")


def flush_position(records, picard_rejected):
    """Group one CHROM/POS bucket by source identity and target END."""
    if not records:
        return

    groups = {}

    for record in records:
        key = merge_key(record)
        groups.setdefault(key, []).append(record)

    # Deterministic output within a genomic position.
    for key in sorted(groups, key=lambda x: (x[2], x[3])):
        flush_group(groups[key], picard_rejected)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--picard-rejected",
        metavar="VCF",
        help="Picard rejected VCF; annotate matching source IDs with rejection metadata",
    )
    args = parser.parse_args()

    picard_rejected = (
        load_picard_rejected(args.picard_rejected)
        if args.picard_rejected
        else {}
    )

    picard_headers_injected = not picard_rejected

    prev_pos = None
    position_records = []

    for raw_line in sys.stdin:
        line = raw_line.rstrip("\n")

        if line.startswith("#"):
            if not picard_headers_injected and line.startswith("#CHROM"):
                for header in _PICARD_HEADERS:
                    sys.stdout.write(header + "\n")
                picard_headers_injected = True

            sys.stdout.write(line + "\n")
            continue

        parts = line.split("\t")

        if len(parts) < 8:
            continue

        pos_key = (parts[0], int(parts[1]))

        if prev_pos is not None and pos_key != prev_pos:
            flush_position(position_records, picard_rejected)
            position_records = []

        position_records.append(parts)
        prev_pos = pos_key

    flush_position(position_records, picard_rejected)


if __name__ == "__main__":
    main()
