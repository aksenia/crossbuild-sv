# crossbuild-sv

A liftover pipeline for structural variants (SVs) that assesses how reliably
hg19 SV calls can be lifted to hg38 — using native hg38 calls on the same
individual as the reference point.

## Motivation

Structural variant calling is reference-dependent. The same individual called
on hg19 and on hg38 natively will produce somewhat different call sets due to
reference differences, caller behaviour, and region accessibility. When only
hg19 calls exist, liftover is the only way to bring them into hg38 coordinates.
But liftover tools disagree, and SVs span large genomic intervals that can
straddle chain-break boundaries — making liftover less reliable than for SNVs.

**This pipeline quantifies that reliability** by comparing liftover results
against native hg38 calls, which — in the absence of a ground truth — serve
as the best available proxy for what the correct hg38 coordinates should be.

## Modules

The project has two modes of operation:

### Per-sample mode

Process one individual to get per-tool concordance metrics against native hg38
calls. This is the core analysis and the output of this pipeline.

```
Same individual
       │
       ├── hg19 SV calls  ──► liftover (4 tools) ──► merged hg38 liftover set
       │                                                          │
       └── hg38 SV calls (native) ────────────────────────────── compare
                                                                  │
                                                            concordance metrics
                                                            per-tool, per-region,
                                                            per-SVTYPE
```

- **hg19 input**: SV VCF from any caller (DRAGEN, Manta, PBSV, …). Filtered
  to target SVTYPE: DEL, DUP, DUP:TANDEM, INS, INV.
- **hg38 input**: ideally the same caller on the same sample natively called on
  hg38, but any hg38 SV call set for the same individual is usable.
- **Liftover**: four independent tools run in parallel on the hg19 VCF. Their
  outputs are merged into an exhaustive set: records where all tools agree are
  collapsed; records where tools disagree on position are kept as separate lines,
  each tagged with which tool produced them.
- **Comparison** (downstream, not in this pipeline): the merged liftover set is
  matched against native hg38 calls by SVTYPE + reciprocal overlap. Concordance
  is evaluated per tool, per genomic region (SEGDUP, CENTEL, CUPs, DISCREPs),
  and per SVTYPE.

### Cohort mode

Run the same per-sample pipeline across multiple individuals, then produce
meta-summaries that aggregate concordance metrics across the cohort. This
answers questions like:

- Which liftover tool is most consistently reliable across samples?
- Are there genomic regions where all tools fail across all samples (systematic
  liftover dead zones)?
- Does per-tool concordance vary by SVTYPE or SV size across the cohort?

```
Sample 1 ──► per-sample pipeline ──► metrics
Sample 2 ──► per-sample pipeline ──► metrics  ──► cohort aggregation ──► meta-summary
Sample N ──► per-sample pipeline ──► metrics
```

Cohort mode is a planned extension. The per-sample pipeline output format is
designed to be directly aggregable: each sample produces a structured metrics
file that the cohort aggregation step can consume without re-running liftover.

## Repository layout

```
sv-preprocess/
├── Dockerfile              # builds crossbuild-sv image (4 liftover tools)
├── liftover/
│   └── merge_tools.py      # merges per-tool outputs on CHROM/POS/ID key
├── snake/
│   ├── Snakefile
│   ├── config.yaml         # user-editable: paths, sv_types, pass_only
│   └── rules/
│       ├── liftover.smk    # filter → stamp IDs → 4 tools → tag → merge
│       └── annotations.smk # region flags on hg19 source + hg38 merged VCF
└── docs/
    ├── pipeline.md         # detailed pipeline steps and how to run
    └── tools.md            # per-tool SV-specific notes and known limitations
```

## Quick start

See [docs/pipeline.md](docs/pipeline.md) for full setup and run instructions.

```bash
# 1. Build container (from repo root)
docker build -f sv-preprocess/Dockerfile -t crossbuild-sv:latest .

# 2. Edit snake/config.yaml with your paths

# 3. Run
snakemake \
    --snakefile snake/Snakefile \
    --configfile snake/config.yaml \
    --use-singularity \
    --singularity-args "..." \
    --cores 8
```

## Outputs

The pipeline produces two annotated VCFs:

| File | Description |
|---|---|
| `annotated/source/<sample>.hg19.annotated.vcf.gz` | Filtered hg19 VCF with region flags; use as input to the comparison tool |
| `annotated/merged/<sample>.merged.hg38.annotated.vcf.gz` | Merged hg38 liftover VCF with `LIFTOVER_TOOL` tags and region flags; compare against native hg38 calls |
