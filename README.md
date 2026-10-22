# crossbuild-sv

A liftover pipeline for structural variants (SVs) that assesses how reliably
hg19 SV calls can be lifted to hg38, using native hg38 calls on the same
individual as the reference point.

## Motivation

Structural variant calling is reference-dependent. The same individual called
on hg19 and on hg38 natively will produce somewhat different call sets due to
reference differences, caller behaviour, and region accessibility. When only
hg19 calls exist, liftover is the only way to bring them into hg38 coordinates.
But liftover tools disagree, and SVs span large genomic intervals that can
straddle chain-break boundaries, making liftover less reliable than for SNVs.

**This pipeline quantifies that reliability** by comparing liftover results
against native hg38 calls, which, in the absence of a ground truth, serve
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
  to target SVTYPEs (default DEL, DUP, INS, INV; DRAGEN tandem duplications
  are `SVTYPE=DUP` with ALT `<DUP:TANDEM>`). BND is out of scope.
- **hg38 input**: ideally the same caller on the same sample natively called on
  hg38, but any hg38 SV call set for the same individual is usable.
- **Liftover**: four independent tools run on the hg19 VCF. Every tool's
  output is normalized the same way (`bcftools norm` + END recomputation), then
  merged into an exhaustive set: records where tools agree on CHROM/POS/END are
  collapsed; disagreements stay as separate lines, each tagged with the tools
  that produced them. Loci where the reference allele changed between builds
  go to a separate `refchanged` set.
- **Comparison** (when `native_vcf` is set): native calls go through the same
  filter and normalization, SVDB links native and lifted calls per tool, and
  per-call links plus long-format metrics (per tool, SVTYPE, size, region,
  precision) are written as flat TSVs. See [docs/comparison.md](docs/comparison.md).

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

Cohort mode is a planned extension. Every per-sample table is long-format
and keyed by `sample_id`, so cohort aggregation is a concatenation (or a
DuckDB glob) of the per-sample files, without re-running liftover.

## Repository layout

```
crossbuild-sv/
├── Dockerfile              # image with the 4 liftover tools, SVDB, truvari, snakemake
├── liftover/
│   ├── finalize_sv.py      # post-norm QC: REF check, END, LIFTOVER_TOOL tag
│   ├── merge_tools.py      # merge per-tool outputs on CHROM/POS/END/ID
│   └── compare_native.py   # links + metrics vs native hg38 calls
├── regions/                # region BEDs (CUPs, DISCREPs, CENTEL, SEGDUP)
├── snake/
│   ├── Snakefile
│   ├── config.yaml         # template: paths, sv_types, contigs, native_vcf
│   └── rules/
│       ├── liftover.smk    # filter → stamp IDs → 4 tools → normalize → merge
│       ├── annotations.smk # region flags on hg19 source + hg38 merged VCFs
│       └── comparison.smk  # native prep → SVDB link → links/metrics
└── docs/
    ├── pipeline.md         # pipeline steps, outputs, how to run
    ├── tools.md            # per-tool SV notes and observed behaviour
    └── comparison.md       # links.tsv.gz / metrics.tsv schema
```

## Quick start

See [docs/pipeline.md](docs/pipeline.md) for full setup and run instructions.

```bash
# 1. Build the container (from the repo root)
docker build -t crossbuild-sv:latest .

# 2. Copy snake/config.yaml and set sample, input_vcf, native_vcf, refs

# 3. Run snakemake inside the container
docker run --rm -v /path/to/ref:/ref:ro -v /path/to/data:/data:ro \
    -v /path/to/results:/results crossbuild-sv:latest \
    snakemake --snakefile /app/snake/Snakefile --configfile /results/config.yaml \
    --directory /results --cores 8
```

## Outputs

| File | Description |
|---|---|
| `annotated/source/<sample>.hg19.annotated.vcf.gz` | Filtered, ID-stamped hg19 VCF with region flags |
| `annotated/merged/<sample>.merged.hg38.annotated.vcf.gz` | Merged hg38 liftover VCF with `LIFTOVER_TOOL` tags and region flags |
| `annotated/merged/<sample>.refchanged.hg38.annotated.vcf.gz` | Same for loci whose reference allele changed between builds (`REFCHANGED_BY`) |
| `comparison/<sample>.links.tsv.gz` | One row per native↔lifted link, per tool and set |
| `comparison/<sample>.metrics.tsv` | Yield and agreement metrics, long format |
