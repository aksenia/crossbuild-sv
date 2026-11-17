#!/usr/bin/env python3
"""Collect and compare allele-state evidence for step 2 of docs/plan.md.

Count unique source IDs, not lines in the ID/tool list.
Keep original genotypes, flags and multiple records visible.
NA means unavailable; none means an inspected record has no such flag.

The classification table is deliberately conservative. Exact allele
comparisons that cannot account for representation changes remain unresolved.
"""

import argparse
import csv
import gzip
from collections import defaultdict
from pathlib import Path


class IndexedFasta:
    """Fetch sequence from an uncompressed FASTA using its .fai index."""

    def __init__(self, path):
        self.path = Path(path)
        self.index = {}
        with open(f"{self.path}.fai") as handle:
            for line in handle:
                fields = line.rstrip("\r\n").split("\t")
                if len(fields) < 5:
                    raise ValueError(f"Malformed FASTA index: {self.path}.fai")
                name, length, offset, line_bases, line_width = fields[:5]
                values = tuple(map(int, (length, offset, line_bases, line_width)))
                if values[0] < 0 or values[1] < 0 or values[2] <= 0 or values[3] < values[2]:
                    raise ValueError(f"Invalid FASTA index entry for {name}: {self.path}.fai")
                self.index[name] = values
        self.handle = open(self.path, "rb")
        if self.handle.read(2) == b"\x1f\x8b":
            self.handle.close()
            raise ValueError("Approach C requires an uncompressed FASTA with a .fai index")

    def close(self):
        self.handle.close()

    def fetch(self, contig, start, end):
        if contig not in self.index:
            raise ValueError(f"Contig {contig!r} is missing from {self.path}.fai")
        length, offset, line_bases, line_width = self.index[contig]
        if start < 0 or end < start or end > length:
            raise ValueError(f"FASTA interval out of range: {contig}:{start}-{end}")
        chunks = []
        cursor = start
        while cursor < end:
            column = cursor % line_bases
            count = min(end - cursor, line_bases - column)
            byte_offset = offset + (cursor // line_bases) * line_width + column
            self.handle.seek(byte_offset)
            chunk = self.handle.read(count)
            if len(chunk) != count:
                raise ValueError(f"Could not read indexed FASTA interval: {contig}:{start}-{end}")
            chunks.append(chunk)
            cursor += count
        return b"".join(chunks).decode("ascii").upper()


def reverse_complement(allele):
    return allele.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1].upper()


def allele_relation(source, lifted_ref):
    """Compare a lifted REF to source alleles without normalizing coordinates."""
    ref = source["ref"].upper()
    alts = source["alt"].upper().split(",")
    lifted_ref = lifted_ref.upper()
    if not lifted_ref or set(lifted_ref) - set("ACGTN"):
        return "unresolved"

    matches_ref = ref in (lifted_ref, reverse_complement(lifted_ref))
    matches_alt = [
        alt in (lifted_ref, reverse_complement(lifted_ref))
        for alt in alts if alt and not (set(alt) - set("ACGTN"))
    ]
    if any(matches_alt) and matches_ref:
        return "ambiguous"
    if any(matches_alt):
        return "swap"
    if matches_ref:
        return "reference_allele"
    if len(matches_alt) == len(alts):
        return "third_allele"
    return "unresolved"


def classification_rows(ids, flagged_by, source, data, fasta):
    rows = []
    for sid in ids:
        tool_flags = []
        swap_reported = False
        raw_relations = []
        fasta_relations = []
        fasta_mismatches = 0

        for tool in ("crossmap", "bcftools", "picard", "transanno"):
            if tool in flagged_by[sid]:
                tool_flags.append(tool)
            raw_records = data[tool]["raw"].get(sid, [])
            for record in raw_records:
                swap = record["SWAP"]
                if swap in ("1", "-1", "present"):
                    swap_reported = True
                    if tool not in tool_flags:
                        tool_flags.append(tool)
                relation = allele_relation(source[sid][0], record["ref"])
                raw_relations.append(relation)

                if set(record["ref"].upper()) - set("ACGTN") or not record["ref"]:
                    continue
                expected = fasta.fetch(
                    record["chrom"], int(record["pos"]) - 1,
                    int(record["pos"]) - 1 + len(record["ref"]),
                )
                if expected != record["ref"].upper():
                    fasta_mismatches += 1
                    continue
                fasta_relations.append(relation)

        approach_a = "swap" if swap_reported else "no_swap_flag"
        decisive_raw = set(raw_relations) - {"unresolved"}
        if "third_allele" in decisive_raw:
            approach_b = "third_allele"
        elif "swap" in decisive_raw and "reference_allele" in decisive_raw:
            approach_b = "conflict"
        elif "swap" in decisive_raw:
            approach_b = "swap"
        elif "ambiguous" in decisive_raw:
            approach_b = "unresolved"
        elif "reference_allele" in decisive_raw:
            approach_b = "reference_allele"
        else:
            approach_b = "unresolved"

        decisive_fasta = set(fasta_relations) - {"unresolved"}
        if (len(decisive_fasta) > 1 or "ambiguous" in decisive_fasta
                or (fasta_mismatches and decisive_fasta)):
            approach_c = "conflict"
        elif decisive_fasta:
            approach_c = next(iter(decisive_fasta))
        elif fasta_mismatches:
            approach_c = "reference_mismatch"
        else:
            approach_c = "unresolved"

        rows.append({
            "source_id": sid,
            "tools": ",".join(tool_flags) or "none",
            "approach_A": approach_a,
            "approach_B": approach_b,
            "approach_C": approach_c,
        })
    return rows


def read_vcf(path, ids, headers):
    records = defaultdict(list)
    opener = gzip.open if path.suffix == ".gz" else open
    samples = []
    with opener(path, "rt") as handle:
        for line in handle:
            if line.startswith(("##INFO=<ID=SWAP,", "##INFO=<ID=REF_CHANGED,")):
                headers.append(f"{path}\n{line.rstrip()}\n")
            if line.startswith("#CHROM"):
                samples = line.rstrip().split("\t")[9:]
                if len(samples) > 1:
                    raise ValueError(f"Expected a single-sample VCF: {path}")
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) < 8:
                raise ValueError(f"Malformed VCF record in {path}")
            if fields[2] not in ids:
                continue
            info = dict(
                (item.split("=", 1) if "=" in item else (item, "present"))
                for item in fields[7].split(";") if item not in ("", ".")
            )
            gt = "NA"
            if samples and len(fields) > 9:
                values = dict(zip(fields[8].split(":"), fields[9].split(":")))
                gt = values.get("GT", "NA")
            records[fields[2]].append({
                "chrom": fields[0], "pos": fields[1],
                "ref": fields[3], "alt": fields[4], "gt": gt,
                "vcf_sample": samples[0] if samples else "NA",
                "filter": fields[6],
                "SWAP": info.get("SWAP", "none"),
                "REF_CHANGED": info.get("REF_CHANGED", "none"),
                "strand_reported": info.get("ORIGINAL_STRAND", "NA"),
                "symbolic_alt": int(any(
                    alt.startswith("<") and alt.endswith(">")
                    for alt in fields[4].split(",")
                )),
                "info": fields[7],
                "extra_columns": "|".join(fields[9 + len(samples):]) or "NA",
            })
    return records


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", type=Path, help="existing pipeline run directory")
    parser.add_argument("--sample", default="HG002", help="pipeline sample name")
    parser.add_argument("--out", type=Path, help="new or empty output directory")
    parser.add_argument("--hg38-fasta", type=Path, required=True,
                        help="uncompressed hg38 FASTA with a .fai index (Approach C)")
    args = parser.parse_args()
    run, sample = args.run.expanduser().resolve(), args.sample
    out = args.out or run / "step2_refstate"

    # 1. find the files
    id_file = run / "liftover/merged" / f"{sample}.refchanged_ids.tsv"
    source_file = run / "liftover/source" / f"{sample}.hg19.vcf.gz"
    tools = ("crossmap", "bcftools", "picard", "transanno")
    paths = {}
    for tool in tools:
        folder = run / "liftover" / tool
        reject = f"{sample}.hg38.unmap.vcf" if tool == "crossmap" else f"{sample}.rejected.vcf"
        paths[tool] = {
            "raw": folder / f"{sample}.hg38.raw.vcf.gz",
            "retained": folder / f"{sample}.hg38.vcf.gz",
            "tool_rejected": folder / reject,
            "qc_rejected": folder / f"{sample}.qc_rejected.vcf",
        }
    inputs = [id_file, source_file]
    inputs.extend(path for stages in paths.values() for path in stages.values())
    for path in inputs:
        if not path.is_file():
            raise ValueError(f"Missing input: {path}")
    if not args.hg38_fasta.is_file() or not Path(f"{args.hg38_fasta}.fai").is_file():
        raise ValueError(f"Approach C requires the FASTA and its .fai index: {args.hg38_fasta}")
    if out.exists() and (not out.is_dir() or any(out.iterdir())):
        raise ValueError(f"Choose a new or empty output directory: {out}")

    # 2. check what we have, using only the refchanged IDs
    flagged_by = defaultdict(set)
    with open(id_file) as handle:
        for line in handle:
            if line.strip():
                sid, tool = line.rstrip("\r\n").split("\t")
                flagged_by[sid].add(tool)
    if not flagged_by:
        raise ValueError("The refchanged ID list is empty")
    ids = sorted(flagged_by)
    headers = []
    source = read_vcf(source_file, flagged_by, headers)
    for sid in ids:
        if len(source[sid]) != 1:
            raise ValueError(f"Expected one source record for {sid}; found {len(source[sid])}")
    data = {
        tool: {stage: read_vcf(path, flagged_by, headers) for stage, path in stages.items()}
        for tool, stages in paths.items()
    }
    summary = [f"Unique refchanged IDs: {len(ids)}", "The plan records 171 for wg4."]
    for label, test in (
        ("Source IDs with symbolic ALT", lambda r: r["symbolic_alt"] == 1),
        ("Source IDs with multiallelic ALT", lambda r: "," in r["alt"]),
        ("Source IDs with sequence-resolved REF and ALTs", lambda r: all(
            allele and set(allele.upper()) <= set("ACGTN")
            for allele in [r["ref"], *r["alt"].split(",")]
        )),
    ):
        summary.append(f"{label}: {sum(test(source[sid][0]) for sid in ids)}")
    for tool in tools:
        counts = {stage: len(records) for stage, records in data[tool].items()}
        missing = sum(not any(data[tool][stage].get(sid) for stage in data[tool]) for sid in ids)
        repeated = sum(len(group) > 1 for records in data[tool].values() for group in records.values())
        summary.append(f"{tool}: {counts}; missing IDs={missing}; repeated ID/stage groups={repeated}")

    # 3. collect the evidence in one table
    # keep a separate row for each record and stage; do not choose a placement
    columns = ["source_id", "tool", "stage", "record_number", "flagged_by",
               "source_chrom", "source_pos", "source_ref", "source_alt", "source_gt",
               "source_vcf_sample", "chrom", "pos", "ref", "alt", "gt", "vcf_sample",
               "filter", "SWAP", "REF_CHANGED", "strand_reported", "symbolic_alt",
               "info", "extra_columns"]
    out.mkdir(parents=True, exist_ok=True)
    with open(out / "evidence.tsv", "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        for sid in ids:
            src = source[sid][0]
            base = {"source_id": sid, "flagged_by": ",".join(sorted(flagged_by[sid]))}
            base.update({f"source_{key}": src[key] for key in ("chrom", "pos", "ref", "alt", "gt", "vcf_sample")})
            for tool in tools:
                found = False
                for stage, records in data[tool].items():
                    for number, record in enumerate(records.get(sid, []), 1):
                        writer.writerow({**base, "tool": tool, "stage": stage,
                                         "record_number": number, **record})
                        found = True
                if not found:
                    row = dict.fromkeys(columns, "NA")
                    writer.writerow({**row, **base, "tool": tool, "stage": "missing"})
    with open(out / "classifications.tsv", "w", newline="") as handle:
        fields = ["source_id", "tools", "approach_A", "approach_B", "approach_C"]
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        fasta = IndexedFasta(args.hg38_fasta)
        try:
            classifications = classification_rows(ids, flagged_by, source, data, fasta)
            writer.writerows(classifications)
        finally:
            fasta.close()
    summary.extend([
        "Counts per stage can overlap; allele categories can overlap.",
        "Approach A: any raw-tool SWAP flag means swap; otherwise no_swap_flag.",
        "Approach B: exact lifted-REF/source-allele comparison; third_allele takes precedence.",
        "Approach C: same comparison after validating lifted REF against the hg38 FASTA.",
        "Exact allele comparisons do not normalize coordinates; unresolved cases need review.",
        f"Classification rows: {len(classifications)} unique source IDs.",
        f"Approach A swaps: {sum(row['approach_A'] == 'swap' for row in classifications)}.",
        f"Approach B swaps: {sum(row['approach_B'] == 'swap' for row in classifications)}.",
        f"Approach C swaps: {sum(row['approach_C'] == 'swap' for row in classifications)}.",
    ])
    (out / "summary.txt").write_text("\n".join(summary) + "\n")
    (out / "flag_headers.txt").write_text("\n".join(headers) or "No flag definitions found.\n")
    print("\n".join(summary))
    print(f"Classifications: {out / 'classifications.tsv'}")
    print(f"Results: {out.resolve()}")


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError) as error:
        raise SystemExit(f"Error: {error}")
