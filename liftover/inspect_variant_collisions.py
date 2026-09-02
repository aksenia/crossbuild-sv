#!/usr/bin/env python3
"""Inspect VCF records that share CHROM/POS/SVTYPE/END.

The script is deliberately read-only.  It answers two separate questions:

1. Which source records have the same coordinate identity key?
2. Are their literal REF/ALT alleles identical, different, or exact reverse
   complements?

For every colliding pair it reports the source VCF line, allele lengths,
first direct difference, common prefix/suffix, reverse-complement tests, the
actual non-matching sequence blocks, QUAL/FILTER, INFO evidence, and every
sample's FORMAT values.  Phasing is assessed only from explicit GT and phase
tags; the script does not infer it from read counts.

Algorithm
---------
1. Read each VCF record and form the key (CHROM, POS, SVTYPE, END).  Missing
   END follows the preprocessing convention and falls back to POS.
2. Retain keys occurring more than once.
3. Compare every pair within each retained key.
4. Test exact equality and exact reverse-complement equality for both the
   complete ALT and the minimal changed ALT segment after common REF/ALT
   prefix/suffix trimming.
5. Use difflib.SequenceMatcher to expose the literal replace/insert/delete
   blocks.  This textual alignment is for inspection, not proof of biological
   equivalence.
6. Report per-record quality, filtering, assembly/read support, genotypes,
   phase sets, and haplotype tags.  Explicit reciprocal 1|0 and 0|1 calls are
   flagged, but only within the same named sample.
"""

from __future__ import annotations

import argparse
import difflib
import gzip
import itertools
import sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path


DNA_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")
DNA_BASES = frozenset("ACGTNacgtn")
SEQUENCE_INFO_FIELDS = frozenset(
    {
        "CONTIG",
        "HOMSEQ",
        "SVINSSEQ",
        "LEFT_SVINSSEQ",
        "RIGHT_SVINSSEQ",
        "DUPSVINSSEQ",
        "DUPHOMSEQ",
    }
)
EVIDENCE_INFO_FIELDS = (
    "CIGAR",
    "CIPOS",
    "CIEND",
    "HOMLEN",
    "HOMSEQ",
    "SVINSLEN",
    "SVINSSEQ",
    "LEFT_SVINSSEQ",
    "RIGHT_SVINSSEQ",
    "DUPSVLEN",
    "DUPSVINSLEN",
    "DUPSVINSSEQ",
    "DUPHOMLEN",
    "DUPHOMSEQ",
    "CONTIG",
    "BND_DEPTH",
    "MATE_BND_DEPTH",
    "JUNCTION_QUAL",
    "EVENT",
    "MATEID",
)
PHASE_FORMAT_FIELDS = ("PS", "HP", "PGT", "PID")


@dataclass(frozen=True)
class SampleCall:
    name: str
    values: dict[str, str]


@dataclass(frozen=True)
class VariantRecord:
    chrom: str
    pos: int
    variant_id: str
    ref: str
    alt: str
    svtype: str
    end: int
    cigar: str
    qual: str
    filters: str
    info: dict[str, str | bool]
    samples: tuple[SampleCall, ...]
    line_number: int
    record_number: int

    @property
    def key(self) -> tuple[str, int, str, int]:
        return (self.chrom, self.pos, self.svtype, self.end)


def open_text(path: str | Path):
    """Open plain VCF or gzip/BGZF-compressed VCF as text."""
    path = str(path)
    if path.endswith((".gz", ".bgz")):
        return gzip.open(path, "rt")
    return open(path, "rt")


def parse_info(text: str) -> dict[str, str | bool]:
    info: dict[str, str | bool] = {}
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


def first_value(value: str | bool | None) -> str | None:
    if value in (None, "", ".", True, False):
        return None
    return str(value).split(",", 1)[0]


def load_records(path: str | Path) -> list[VariantRecord]:
    records: list[VariantRecord] = []
    record_number = 0
    sample_names: list[str] = []

    with open_text(path) as handle:
        for line_number, raw in enumerate(handle, start=1):
            if raw.startswith("#CHROM"):
                header_parts = raw.rstrip("\n").split("\t")
                sample_names = header_parts[9:]
                continue
            if raw.startswith("#") or not raw.strip():
                continue

            record_number += 1
            parts = raw.rstrip("\n").split("\t")
            if len(parts) < 8:
                raise ValueError(
                    f"Malformed VCF record at decompressed line {line_number}: "
                    f"expected at least 8 columns"
                )

            info = parse_info(parts[7])
            pos = int(parts[1])
            end_text = first_value(info.get("END"))
            end = int(end_text) if end_text is not None else pos
            svtype = first_value(info.get("SVTYPE")) or "."
            cigar = first_value(info.get("CIGAR")) or "."
            format_keys = (
                parts[8].split(":")
                if len(parts) > 8 and parts[8] not in ("", ".")
                else []
            )
            sample_calls: list[SampleCall] = []
            for sample_index, sample_text in enumerate(parts[9:]):
                sample_name = (
                    sample_names[sample_index]
                    if sample_index < len(sample_names)
                    else f"SAMPLE_{sample_index + 1}"
                )
                sample_values = sample_text.split(":")
                values = {
                    key: sample_values[index] if index < len(sample_values) else "."
                    for index, key in enumerate(format_keys)
                }
                sample_calls.append(SampleCall(sample_name, values))

            records.append(
                VariantRecord(
                    chrom=parts[0],
                    pos=pos,
                    variant_id=parts[2],
                    ref=parts[3],
                    alt=parts[4],
                    svtype=svtype,
                    end=end,
                    cigar=cigar,
                    qual=parts[5],
                    filters=parts[6],
                    info=info,
                    samples=tuple(sample_calls),
                    line_number=line_number,
                    record_number=record_number,
                )
            )

    return records


def collision_groups(
    records: list[VariantRecord],
    positions: set[tuple[str, int]] | None = None,
) -> list[tuple[tuple[str, int, str, int], list[VariantRecord]]]:
    grouped: dict[tuple[str, int, str, int], list[VariantRecord]] = defaultdict(list)
    for record in records:
        if positions and (record.chrom, record.pos) not in positions:
            continue
        grouped[record.key].append(record)

    return sorted(
        ((key, values) for key, values in grouped.items() if len(values) > 1),
        key=lambda item: (item[0][0], item[0][1], item[0][2], item[0][3]),
    )


def is_literal_dna(sequence: str) -> bool:
    return bool(sequence) and all(base in DNA_BASES for base in sequence)


def reverse_complement(sequence: str) -> str | None:
    if not is_literal_dna(sequence):
        return None
    return sequence.translate(DNA_COMPLEMENT)[::-1].upper()


def common_prefix_length(first: str, second: str) -> int:
    limit = min(len(first), len(second))
    index = 0
    while index < limit and first[index] == second[index]:
        index += 1
    return index


def common_suffix_length(first: str, second: str, prefix_length: int = 0) -> int:
    limit = min(len(first), len(second)) - prefix_length
    index = 0
    while index < limit and first[-index - 1] == second[-index - 1]:
        index += 1
    return index


def minimal_changed_alt(ref: str, alt: str) -> str | None:
    """Return ALT after removing shared REF/ALT prefix and suffix.

    This exposes the sequence introduced by the allele without assuming that
    INFO/CIGAR is present.  It is a literal trimming operation, not a genomic
    normalization algorithm.
    """
    if not (is_literal_dna(ref) and is_literal_dna(alt)):
        return None

    prefix = common_prefix_length(ref, alt)
    suffix = common_suffix_length(ref, alt, prefix)
    stop = len(alt) - suffix if suffix else len(alt)
    return alt[prefix:stop].upper()


def exact_reverse_complement(first: str | None, second: str | None) -> bool | None:
    if first is None or second is None:
        return None
    second_rc = reverse_complement(second)
    if second_rc is None:
        return None
    return first.upper() == second_rc


def direct_hamming(first: str, second: str) -> int | None:
    if len(first) != len(second):
        return None
    return sum(left != right for left, right in zip(first, second))


def format_boolean(value: bool | None) -> str:
    if value is None:
        return "not-applicable"
    return "YES" if value else "no"


def sequence_preview(sequence: str, flank: int, show_full: bool) -> str:
    if sequence == "":
        return "<empty>"
    if show_full or len(sequence) <= 2 * flank + 12:
        return sequence
    omitted = len(sequence) - 2 * flank
    return f"{sequence[:flank]}...<{omitted} bp omitted>...{sequence[-flank:]}"


def present_sample_value(value: str | None) -> bool:
    return value not in (None, "", ".")


def genotype_phase(gt: str | None) -> str:
    if not present_sample_value(gt):
        return "not available"
    if "|" in str(gt):
        return "phased"
    if "/" in str(gt):
        return "unphased"
    return "phase delimiter absent"


def reciprocal_phased_genotypes(first_gt: str | None, second_gt: str | None) -> bool | None:
    if not (present_sample_value(first_gt) and present_sample_value(second_gt)):
        return None
    return (str(first_gt), str(second_gt)) in {
        ("1|0", "0|1"),
        ("0|1", "1|0"),
    }


def sample_map(record: VariantRecord) -> dict[str, SampleCall]:
    return {sample.name: sample for sample in record.samples}


def format_info_evidence_value(
    field: str,
    value: str | bool,
    context: int,
    show_full: bool,
) -> str:
    if value is True:
        return "present"
    text = str(value)
    if field in SEQUENCE_INFO_FIELDS:
        return f"{len(text)} bp; {sequence_preview(text, context, show_full)}"
    return text


def print_record_evidence(
    label: str,
    record: VariantRecord,
    context: int,
    show_full_evidence: bool,
) -> None:
    print(f"    {label} evidence: QUAL={record.qual}; FILTER={record.filters}")

    found_info = False
    for field in EVIDENCE_INFO_FIELDS:
        value = record.info.get(field)
        if value in (None, "", ".", False):
            continue
        found_info = True
        formatted = format_info_evidence_value(
            field, value, context, show_full_evidence
        )
        print(f"      INFO/{field}={formatted}")
    if not found_info:
        print("      INFO evidence fields: none present")

    if not record.samples:
        print("      Sample FORMAT evidence: no sample columns")
        return

    for sample in record.samples:
        gt = sample.values.get("GT")
        format_text = "; ".join(
            f"{field}={value}" for field, value in sample.values.items()
        )
        print(
            f"      Sample {sample.name}: GT={gt or '.'} "
            f"({genotype_phase(gt)})"
        )
        print(f"        FORMAT: {format_text or 'none'}")
        phase_tags = [
            f"{field}={sample.values[field]}"
            for field in PHASE_FORMAT_FIELDS
            if present_sample_value(sample.values.get(field))
        ]
        print(
            "        Explicit phase/haplotype tags: "
            f"{'; '.join(phase_tags) if phase_tags else 'none present'}"
        )


def print_pair_phasing(first: VariantRecord, second: VariantRecord) -> None:
    first_samples = sample_map(first)
    second_samples = sample_map(second)
    common_samples = sorted(set(first_samples).intersection(second_samples))

    print("    Pair phasing assessment:")
    if not common_samples:
        print("      No common named samples with FORMAT values")
        return

    for sample_name in common_samples:
        first_call = first_samples[sample_name]
        second_call = second_samples[sample_name]
        first_gt = first_call.values.get("GT")
        second_gt = second_call.values.get("GT")
        reciprocal = reciprocal_phased_genotypes(first_gt, second_gt)
        first_ps = first_call.values.get("PS")
        second_ps = second_call.values.get("PS")

        print(
            f"      {sample_name}: A GT={first_gt or '.'} "
            f"({genotype_phase(first_gt)}); B GT={second_gt or '.'} "
            f"({genotype_phase(second_gt)})"
        )
        print(
            "        Reciprocal 1|0 versus 0|1: "
            f"{format_boolean(reciprocal)}"
        )
        if present_sample_value(first_ps) and present_sample_value(second_ps):
            relation = "same" if first_ps == second_ps else "different"
            print(
                f"        PS: A={first_ps}; B={second_ps}; relation={relation}"
            )
        else:
            print(
                f"        PS: A={first_ps or 'not available'}; "
                f"B={second_ps or 'not available'}"
            )


def print_pair_sequence_evidence(first: VariantRecord, second: VariantRecord) -> None:
    comparisons = []
    for field in EVIDENCE_INFO_FIELDS:
        if field not in SEQUENCE_INFO_FIELDS:
            continue
        first_value_raw = first.info.get(field)
        second_value_raw = second.info.get(field)
        if not isinstance(first_value_raw, str) or not isinstance(second_value_raw, str):
            continue
        if first_value_raw in ("", ".") or second_value_raw in ("", "."):
            continue
        comparisons.append((field, first_value_raw, second_value_raw))

    print("    Pair assembly/sequence evidence comparison:")
    if not comparisons:
        print("      No shared sequence-valued INFO evidence fields")
        return

    for field, first_sequence, second_sequence in comparisons:
        prefix = common_prefix_length(first_sequence, second_sequence)
        first_difference = None if first_sequence == second_sequence else prefix + 1
        hamming = direct_hamming(first_sequence, second_sequence)
        print(
            f"      INFO/{field}: lengths={len(first_sequence)}/{len(second_sequence)}; "
            f"exact={format_boolean(first_sequence == second_sequence)}; "
            "reverse-complement="
            f"{format_boolean(exact_reverse_complement(first_sequence, second_sequence))}; "
            "hamming="
            f"{hamming if hamming is not None else 'not-applicable'}; "
            f"first difference={first_difference}"
        )


def edit_blocks(first: str, second: str) -> list[tuple[str, int, int, int, int]]:
    matcher = difflib.SequenceMatcher(None, first, second, autojunk=False)
    return [block for block in matcher.get_opcodes() if block[0] != "equal"]


def format_interval(start: int, stop: int) -> str:
    if start == stop:
        return f"empty before base {start + 1}"
    return f"bases {start + 1}-{stop}"


def print_edit_blocks(
    first: str,
    second: str,
    context: int,
    max_blocks: int,
    show_full: bool,
) -> None:
    blocks = edit_blocks(first, second)
    visible = blocks if max_blocks == 0 else blocks[:max_blocks]

    if not blocks:
        print("    Sequence change blocks: none")
        return

    print(f"    Sequence change blocks: {len(blocks)}")
    for number, (operation, first_start, first_stop, second_start, second_stop) in enumerate(
        visible, start=1
    ):
        first_change = first[first_start:first_stop]
        second_change = second[second_start:second_stop]
        print(f"      {number}. {operation.upper()}")
        print(
            f"         A {format_interval(first_start, first_stop)} "
            f"({len(first_change)} bp): "
            f"{sequence_preview(first_change, context, show_full)}"
        )
        print(
            f"         B {format_interval(second_start, second_stop)} "
            f"({len(second_change)} bp): "
            f"{sequence_preview(second_change, context, show_full)}"
        )

    if len(visible) < len(blocks):
        print(
            f"      ... {len(blocks) - len(visible)} additional blocks hidden; "
            "use --max-blocks 0 to show all"
        )


def print_pair(
    first: VariantRecord,
    second: VariantRecord,
    context: int,
    max_blocks: int,
    show_full: bool,
    show_full_evidence: bool,
) -> None:
    alt_prefix = common_prefix_length(first.alt, second.alt)
    alt_suffix = common_suffix_length(first.alt, second.alt, alt_prefix)
    first_difference = None if first.alt == second.alt else alt_prefix + 1
    first_core = minimal_changed_alt(first.ref, first.alt)
    second_core = minimal_changed_alt(second.ref, second.alt)
    hamming = direct_hamming(first.alt, second.alt)

    print(f"  Pair: {first.variant_id}  <->  {second.variant_id}")
    print(
        f"    A: data record {first.record_number}, decompressed line "
        f"{first.line_number}, CIGAR={first.cigar}"
    )
    print(
        f"    B: data record {second.record_number}, decompressed line "
        f"{second.line_number}, CIGAR={second.cigar}"
    )
    print(f"    REF: {first.ref}  <->  {second.ref}")
    print(f"    ALT lengths: {len(first.alt)}  <->  {len(second.alt)}")
    print(f"    Exact REF equality: {format_boolean(first.ref == second.ref)}")
    print(f"    Exact ALT equality: {format_boolean(first.alt == second.alt)}")
    print(
        "    ALT A == reverse-complement(ALT B): "
        f"{format_boolean(exact_reverse_complement(first.alt, second.alt))}"
    )
    print(
        "    Changed ALT A == reverse-complement(changed ALT B): "
        f"{format_boolean(exact_reverse_complement(first_core, second_core))}"
    )
    print(
        "    Direct Hamming distance: "
        f"{hamming if hamming is not None else 'not-applicable (different lengths)'}"
    )
    print(f"    Common ALT prefix: {alt_prefix} bp")
    print(f"    Common ALT suffix: {alt_suffix} bp")
    print(f"    First direct ALT difference: {first_difference}")
    print(
        "    Minimal changed ALT lengths: "
        f"{len(first_core) if first_core is not None else 'not-applicable'}  <->  "
        f"{len(second_core) if second_core is not None else 'not-applicable'}"
    )

    print_record_evidence(
        "A", first, context=context, show_full_evidence=show_full_evidence
    )
    print_record_evidence(
        "B", second, context=context, show_full_evidence=show_full_evidence
    )
    print_pair_sequence_evidence(first, second)
    print_pair_phasing(first, second)

    if show_full:
        print(f"    Full A REF: {first.ref}")
        print(f"    Full A ALT: {first.alt}")
        print(f"    Full B REF: {second.ref}")
        print(f"    Full B ALT: {second.alt}")

    print_edit_blocks(first.alt, second.alt, context, max_blocks, show_full)


def parse_position(text: str) -> tuple[str, int]:
    try:
        chrom, position = text.rsplit(":", 1)
        return chrom, int(position.replace(",", ""))
    except (ValueError, TypeError) as error:
        raise argparse.ArgumentTypeError(
            f"expected CHROM:POS, received {text!r}"
        ) from error


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Find VCF records sharing CHROM/POS/SVTYPE/END and show their "
            "literal allele and reverse-complement differences."
        )
    )
    parser.add_argument("vcf", help="Source .vcf, .vcf.gz, or .vcf.bgz file")
    parser.add_argument(
        "--position",
        action="append",
        type=parse_position,
        metavar="CHROM:POS",
        help="inspect only this position; may be supplied more than once",
    )
    parser.add_argument(
        "--context",
        type=int,
        default=25,
        help="bases shown at each end of a long change block (default: 25)",
    )
    parser.add_argument(
        "--max-blocks",
        type=int,
        default=20,
        help="maximum change blocks per pair; 0 shows all (default: 20)",
    )
    parser.add_argument(
        "--show-full-alleles",
        action="store_true",
        help="print complete REF and ALT strings in addition to change blocks",
    )
    parser.add_argument(
        "--show-full-evidence-sequences",
        action="store_true",
        help="print complete sequence-valued INFO evidence such as CONTIG",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.context < 1:
        raise SystemExit("--context must be at least 1")
    if args.max_blocks < 0:
        raise SystemExit("--max-blocks must be 0 or greater")

    records = load_records(args.vcf)
    positions = set(args.position) if args.position else None
    groups = collision_groups(records, positions)

    pair_count = sum(len(group) * (len(group) - 1) // 2 for _, group in groups)
    record_count = sum(len(group) for _, group in groups)
    print(f"VCF: {args.vcf}")
    print("Collision key: CHROM + POS + SVTYPE + END")
    print(
        f"Found {len(groups)} collision group(s), {record_count} record(s), "
        f"and {pair_count} pairwise comparison(s)."
    )

    for group_number, (key, group) in enumerate(groups, start=1):
        chrom, pos, svtype, end = key
        print()
        print(
            f"=== Group {group_number}: CHROM={chrom} POS={pos} "
            f"SVTYPE={svtype} END={end} ({len(group)} records) ==="
        )
        for first, second in itertools.combinations(group, 2):
            print_pair(
                first,
                second,
                context=args.context,
                max_blocks=args.max_blocks,
                show_full=args.show_full_alleles,
                show_full_evidence=args.show_full_evidence_sequences,
            )

    if not groups:
        print("No duplicate coordinate keys were found for the selected records.")

    return 0


if __name__ == "__main__":
    sys.exit(main())
