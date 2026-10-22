# Pipeline

Lifts a single-sample hg19 SV VCF to hg38 with four tools, normalizes and
merges their outputs into one exhaustive hg38 set, adds region flags, and
(optionally) compares the lifted calls with native hg38 calls of the same
individual.

## Overview

```
input hg19 SV VCF
       │
       ▼
 prep_source_vcf        filter SVTYPE (sv_types), PASS (pass_only), contigs;
                        stamp ID = CHROM_POS_SVTYPE_END (hg19), _2,_3 on collisions
       │
       ├──────────────────────────────┐
       ▼                              ▼
 normalize_contigs              transanno liftvcf
 (add chr prefix)               (bare contigs, own chain)
       │
       ├────────────┬────────────────┐
       ▼            ▼                ▼
 CrossMap vcf  bcftools +liftover  Picard LiftoverVcf
               (--no-left-align)   (LIFTOVER_MIN_MATCH=0.0)
       │            │                │               │
       └────────────┴───────┬────────┴───────────────┘
                            ▼
 normalize_liftover     per tool: bcftools norm (non-symbolic records) +
                        finalize_sv.py (REF check, END, LIFTOVER_TOOL) → qc_rejected
                            │
                            ▼
 refchanged_ids         IDs any tool flags as REF-changed (multiallelic / SWAP / REF_CHANGED)
                            │
                            ▼
 merge_liftovers        collapse identical CHROM/POS/END/ID; disagreements stay
                        as separate lines; REF-changed IDs → refchanged set
                            │
               ┌────────────┴─────────────┐
               ▼                          ▼
 annotate_regions_hg19         annotate_regions_hg38
 (CUPs, DISCREPs, CENTEL,      (CUPs, DISCREPs; merged + refchanged)
  SEGDUP)
                            │
                            ▼  (only when native_vcf is set)
 comparison             native prep → per tool × set SVDB link → links + metrics
                        (docs/comparison.md)
```

### Key design choices

- **ID scheme**: each source record is stamped `CHROM_POS_SVTYPE_END`
  (hg19 coords) with awk before liftover (`bcftools --set-id` misparses
  `%INFO/SVTYPE_`). Distinct records sharing a key get `_2`, `_3`, … in file
  order. All tools preserve the ID, so a source variant is tracked across
  tools regardless of its hg38 coordinates. Native hg38 calls get the same
  stamp in hg38 coordinates.
- **Normalization before merge**: tools differ in allele spelling (bcftools
  pads REF/ALT, transanno drops a shared padding base). `bcftools norm`
  left-aligns and trims sequence-resolved records; `finalize_sv.py` checks REF
  against hg38 and sets END from the normalized record. Symbolic records
  bypass `bcftools norm` (no sequence to align; norm 1.22 shifts
  `<DUP:TANDEM>` and resets END=POS); their END is mapped through the chain
  from the source POS/END, accepted only if the chain maps the source POS onto
  the tool's POS.
- **Merge key**: `CHROM/POS/END/ID` (not REF/ALT; SV alleles can be kb
  long). Tools agreeing on a placement collapse into one record with
  comma-joined `LIFTOVER_TOOL`; any remaining difference is a real
  disagreement and stays on its own line.
- **REF-changed set**: at loci where the hg38 reference carries a different
  allele than hg19 (multiallelic ALT, bcftools `SWAP`, transanno
  `REF_CHANGED`) the lifted record is not a like-for-like SV. Every tool's
  record of such an ID goes to `merged/<s>.refchanged.hg38.vcf.gz` with
  `REFCHANGED_BY`, so a tool that silently mis-handles the locus does not
  leave it in the main set.

## Supported SV types

Configured via `sv_types`; anything else is dropped at the filter step. BND is
out of scope.

| SVTYPE | Notes |
|---|---|
| DEL | Sequence-resolved or symbolic `<DEL>` |
| DUP | Symbolic `<DUP>`; DRAGEN tandem duplications are `SVTYPE=DUP`, ALT `<DUP:TANDEM>` |
| INS | Sequence-resolved or symbolic `<INS>` |
| INV | Symbolic `<INV>` |

## Setup

### 1. Build the container

```bash
docker build -t crossbuild-sv:latest .
```

The image contains CrossMap 0.7.3, bcftools 1.22 + liftover plugin, Picard
3.4.0, transanno 0.4.5, SVDB 2.12.0, truvari 5.3.0, snakemake 8.30, the
pipeline code under `/app` and the region BEDs under `/app/regions/`.

For HPC / Apptainer:

```bash
apptainer build crossbuild-sv.sif docker-daemon://crossbuild-sv:latest
```

### 2. Reference data

- **hg19 FASTA** + `.fai`, numeric contigs (`1`, `2`, …): Picard source, transanno
- **hg19 FASTA chr-prefixed** + `.fai`: `bcftools +liftover` source
- **hg38 FASTA** + `.fai` + `.dict`, chr-prefixed
- **hg19→hg38 UCSC chain** (`hg19ToHg38.over.chain`): CrossMap, bcftools, Picard
- **Transanno chain**: built with minimap2 (`minimap2 -cx asm5 --cs hg38.fa
  hg19.fa > paf; transanno minimap2chain paf --output chain`), contigs named
  like the numeric hg19 FASTA

### 3. Configure

Copy `snake/config.yaml` and set at least `sample`, `input_vcf`,
`results_dir` and the `ref:` paths (all as seen inside the container). There
is no default configfile: `--configfile` is required.

| Key | Meaning |
|---|---|
| `sample` | sample ID, used in file names and `sample_id` columns |
| `input_vcf` | hg19 SV VCF (single sample) |
| `native_vcf` | optional native hg38 VCF of the same individual; enables comparison |
| `contigs` | optional hg19 contig subset in source naming (e.g. `["14", "21", "Y"]`); `[]` = whole genome. Native is restricted to the same contigs |
| `sv_types`, `pass_only` | source/native filter |
| `liftover_tools` | subset/order of `crossmap, bcftools, picard, transanno` |
| `comparison.svdb_overlap`, `comparison.svdb_profile` | SVDB linking parameters |

### 4. Run

Snakemake runs inside the container (all tools are there):

```bash
docker run --rm --user $(id -u):$(id -g) --env HOME=/tmp \
    -v /path/to/ref:/ref:ro -v /path/to/data:/data:ro \
    -v /path/to/results:/results crossbuild-sv:latest \
    snakemake --snakefile /app/snake/Snakefile --configfile /results/config.yaml \
    --directory /results --cores 8
```

Add `-n` for a dry run. If reference files are symlinks, mount their targets
at the same path. To work on the code without rebuilding, mount the
repo's `snake/`, `liftover/` and `regions/` over `/app/...`.

With Apptainer, the same command works via `apptainer exec -B ... crossbuild-sv.sif snakemake ...`.

A whole-genome HG002 DRAGEN run (~17k source SVs, with comparison) takes
about 2–3 minutes on 8 cores.

## Outputs

| File | Description |
|---|---|
| `liftover/source/<s>.hg19.vcf.gz` | Filtered, ID-stamped source VCF |
| `liftover/<tool>/<s>.hg38.raw.vcf.gz` | Tool output as written by the tool (sorted) |
| `liftover/<tool>/<s>.hg38.vcf.gz` | Normalized, finalized per-tool hg38 VCF (`LIFTOVER_TOOL`) |
| `liftover/<tool>/<s>.rejected.vcf` | Tool rejections (CrossMap: `<s>.hg38.unmap.vcf`) |
| `liftover/<tool>/<s>.qc_rejected.vcf` | Lifted but failed post-liftover QC (`QC_REJECT_REASON`) |
| `liftover/<tool>/<s>.normalize.log` | `bcftools norm` / finalize log |
| `liftover/transanno/<s>.log` | Transanno run log |
| `liftover/merged/<s>.refchanged_ids.tsv` | `ID<TAB>tool` for every REF-changed flag |
| `annotated/source/<s>.hg19.annotated.vcf.gz` | Source VCF with region flags |
| `annotated/merged/<s>.merged.hg38.annotated.vcf.gz` | Merged main set with region flags |
| `annotated/merged/<s>.refchanged.hg38.annotated.vcf.gz` | Merged REF-changed set with region flags |
| `comparison/native/<s>.hg38.vcf.gz` | Prepared native calls (+ `qc_rejected.vcf`) |
| `comparison/svdb/<tool>.<set>.vcf` | Raw SVDB link output |
| `comparison/<s>.links.tsv.gz`, `comparison/<s>.metrics.tsv` | See [comparison.md](comparison.md) |

### INFO fields in the merged hg38 VCFs

| Field | Description |
|---|---|
| `LIFTOVER_TOOL` | Comma-joined tools that produced this placement |
| `END` | hg38 end, recomputed by `finalize_sv.py` |
| `SVTYPE`, `SVLEN` | From the source record |
| `PICARD_REJECT_REASON` (+ `PICARD_ATTEMPTED_*`) | Why Picard rejected this source variant |
| `MULTIMAP`, `REF_CHANGED` | Carried from transanno |
| `REFCHANGED_BY` | refchanged set: tools that flagged the ID |
| `CUPs` / `DISCREPs` | hg38 region overlap flags (`CENTEL`, `SEGDUP` are on the hg19 source only) |
