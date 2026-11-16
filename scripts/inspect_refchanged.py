#!/usr/bin/env python3
"""Collect evidence for step 2 of docs/plan.md. No classification yet.

Count unique source IDs, not lines in the ID/tool list.
Keep original genotypes, flags and multiple records visible.
NA means unavailable; none means an inspected record has no such flag.
"""

import argparse
import csv
import gzip
from collections import defaultdict
from pathlib import Path


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
    summary.extend(["Counts per stage can overlap; allele categories can overlap.",
                    "No classification, normalization or hg38 REF check performed."])
    (out / "summary.txt").write_text("\n".join(summary) + "\n")
    (out / "flag_headers.txt").write_text("\n".join(headers) or "No flag definitions found.\n")
    print("\n".join(summary))
    print(f"Results: {out.resolve()}")


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError) as error:
        raise SystemExit(f"Error: {error}")
