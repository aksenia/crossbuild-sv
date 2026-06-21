# sv-preprocess

Liftover pipeline for structural variants (SVs) from hg19 to hg38. Runs four
independent tools in parallel, merges their outputs, and annotates with
genomic region flags (CUPs, DISCREPs, CENTEL, SEGDUP).

Designed as a companion to `preprocess/` (which handles SNVs/indels). The
output VCFs feed into a downstream SV comparison and prioritisation tool.

## Pipeline overview

```
input hg19 SV VCF
       │
       ▼
 filter to target SVTYPE           config: sv_types (DEL, DUP, DUP:TANDEM, INS, INV)
 optionally require PASS            config: pass_only (default: true)
 stamp hg19 ID: CHROM_POS_SVTYPE_END
       │
       ├──────────────────────────────────────────────┐
       │                                              │
       ▼                                              ▼
 normalize contigs                           transanno liftvcf
 (add chr prefix)                            (uses bare contig names)
       │
       ├──────────────┬─────────────────┐
       ▼              ▼                 ▼
 CrossMap vcf   bcftools +liftover   Picard LiftoverVcf
                (--no-left-align)    (LIFTOVER_MIN_MATCH=0.0)
       │              │                 │               │ (transanno)
       └──────────────┴─────────────────┴───────────────┘
                             │
                             ▼
                  tag each output with LIFTOVER_TOOL
                             │
                             ▼
                   merge (collapse identical
                   CHROM/POS/ID records;
                   keep separate lines when
                   tools disagree on position)
                             │
               ┌─────────────┴─────────────┐
               ▼                           ▼
  annotate hg19 source VCF      annotate merged hg38 VCF
  (CUPs, DISCREPs, CENTEL,      (CUPs, DISCREPs)
   SEGDUP)
```

### Key design choices

- **ID scheme**: each variant is stamped `CHROM_POS_SVTYPE_END` (hg19 coords)
  before liftover. All tools preserve the ID, so concordance can be tracked
  across tools regardless of hg38 coordinate.
- **Merge key**: records are merged on `CHROM/POS/ID` (hg38). Tools that lift
  the same variant to the same position are collapsed into one record with a
  comma-joined `LIFTOVER_TOOL` INFO value. Tools that disagree on position
  produce separate records.
- **No VEP**: VEP annotation is out of scope for this pipeline.

## Supported SV types

Configured via `sv_types` in `config.yaml`. Anything else (e.g. BND) is
dropped at the filter step.

| SVTYPE | Notes |
|---|---|
| DEL | Sequence-resolved or symbolic `<DEL>` |
| DUP | Symbolic `<DUP>` |
| DUP:TANDEM | DRAGEN tandem duplication |
| INS | Sequence-resolved or symbolic `<INS>` |
| INV | Symbolic `<INV>` |

## Setup

### 1. Build the container

Build from the **repo root** (so the Dockerfile can reference `preprocess/regions/`):

```bash
docker build -f sv-preprocess/Dockerfile -t crossbuild-sv:latest .
```

For HPC / air-gapped environments, export as a Singularity/Apptainer SIF:

```bash
apptainer build crossbuild-sv.sif docker-daemon://crossbuild-sv:latest
```

### 2. Extract pipeline files from the SIF (one-time, on HPC)

Snakemake must run **on the host** so it can dispatch per-rule containers.
Extract `snake/` from the SIF once:

```bash
cd /your/analysis/dir
singularity exec crossbuild-sv.sif bash -c 'cd /app && tar cf - snake' \
    | tar -xf - -C .
```

### 3. Configure

Edit `snake/config.yaml` (extracted from the SIF or copied from this repo):

```yaml
sample: my_sample
input_vcf: /data/my_sample.hg19.sv.vcf.gz
results_dir: /results

sv_types: [DEL, DUP, "DUP:TANDEM", INS, INV]
pass_only: true          # set false to keep non-PASS records

liftover_tools: [crossmap, bcftools, picard, transanno]

ref:
  hg19_fasta: /ref/hg19.fa
  hg19_fasta_chr: /ref/hg19.chr.fa     # chr-prefixed, required by bcftools +liftover
  hg38_fasta: /ref/hg38.fa
  liftover_chain: /ref/hg19ToHg38.over.chain
  transanno_chain: /ref/transanno.hg19tohg38.chain

container_image: "/path/to/crossbuild-sv.sif"   # or docker://crossbuild-sv:latest
```

Region BED paths default to `/app/regions/` inside the container — no changes
needed unless you bind-mount them elsewhere.

### 4. Run

```bash
module load Singularity
module load snakemake/8.4.2

snakemake \
    --snakefile /path/to/snake/Snakefile \
    --configfile /path/to/config.yaml \
    --use-singularity \
    --singularity-args "--no-home \
        -B /path/to/hg19.fa:/ref/hg19.fa \
        -B /path/to/hg19.fa.fai:/ref/hg19.fa.fai \
        -B /path/to/hg19.chr.fa:/ref/hg19.chr.fa \
        -B /path/to/hg19.chr.fa.fai:/ref/hg19.chr.fa.fai \
        -B /path/to/hg38.fa:/ref/hg38.fa \
        -B /path/to/hg38.fa.fai:/ref/hg38.fa.fai \
        -B /path/to/hg38.fa.dict:/ref/hg38.fa.dict \
        -B /path/to/hg19ToHg38.over.chain:/ref/hg19ToHg38.over.chain \
        -B /path/to/transanno.chain:/ref/transanno.hg19tohg38.chain \
        -B /path/to/data:/data \
        -B /path/to/results:/results" \
    --cores 8
```

Add `-n` for a dry run.

## Reference data

Same requirements as `preprocess/` for the liftover references:

- **hg19 FASTA** + `.fai` — numeric contig names (`1`, `2` …)
- **hg19 FASTA chr-prefixed** + `.fai` — required by `bcftools +liftover`
- **hg38 FASTA** + `.fai` + `.dict`
- **hg19→hg38 UCSC chain** (`hg19ToHg38.over.chain`)
- **Transanno chain** — minimap2-generated; see `preprocess/` docs for build instructions

All region BED files (CUPs, DISCREPs, CENTEL, SEGDUP) are pre-built and
shipped inside the container at `/app/regions/`.

## Outputs

| File | Description |
|---|---|
| `annotated/source/<sample>.hg19.annotated.vcf.gz` | Filtered, ID-stamped hg19 VCF with region flags |
| `annotated/merged/<sample>.merged.hg38.annotated.vcf.gz` | Merged hg38 VCF with LIFTOVER_TOOL + region flags |
| `liftover/<tool>/<sample>.hg38.vcf.gz` | Per-tool hg38 liftover output |
| `liftover/<tool>/<sample>.rejected.vcf` | Per-tool rejected records |
| `liftover/transanno/<sample>.log` | Transanno run log |

### Key INFO fields in the merged hg38 VCF

| Field | Description |
|---|---|
| `LIFTOVER_TOOL` | Comma-joined list of tools that produced this position (e.g. `crossmap,bcftools,picard`) |
| `SVTYPE` | Original SV type from hg19 |
| `END` | hg38 end coordinate (updated by all tools except where noted) |
| `SVLEN` | SV length (preserved from hg19) |
| `PICARD_REJECT_REASON` | Why Picard rejected this variant (when applicable) |
| `CUPs` / `DISCREPs` / `CENTEL` / `SEGDUP` | Region overlap flags (hg19 source VCF) |
