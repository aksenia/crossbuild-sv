# Liftover tools used by crossbuild-sv

The pipeline runs four independent tools on the same filtered, ID-stamped
source records. This document lists the commands that are actually part of the
workflow. Optional liftover modes are intentionally not added.

## CrossMap

Container command:

```text
CrossMap vcf --chromid l <UCSC-chain> <chr-hg19.vcf.gz> <hg38.fa> <raw.vcf>
```

`--chromid l` keeps the long, chr-prefixed target names used by hg38, the
target BED files, and the other branches. CrossMap writes unmapped records to
its `.unmap` output, which the rule retains as
`<sample>.hg38.unmap.vcf`.

The workflow does not use `CrossMap region` and does not apply a separate
chain-coverage filter.

## bcftools +liftover

Container version: bcftools 1.22 with the freeseek liftover plugin.

```text
bcftools +liftover -Ov <chr-hg19.vcf.gz> -- \
  -c <UCSC-chain> \
  --src-fasta-ref <chr-hg19.fa> \
  --fasta-ref <hg38.fa> \
  --no-left-align \
  --reject <rejected.vcf> \
  --write-reject
```

`--no-left-align` prevents this branch from adding a post-liftover
left-alignment step. The rule does not use `--lift-end`; the same explicit
END-only contract is applied after every tool.

## Picard LiftoverVcf

Container version: Picard 3.4.0.

Before execution, only the VCF fileformat declaration is reduced to VCFv4.3
because Picard does not accept a newer declaration. Records and alleles are not
pre-filtered.

```text
JAVA_TOOL_OPTIONS=-Xmx6g picard LiftoverVcf \
  -INPUT <chr-hg19.vcf.gz> \
  -OUTPUT <raw.vcf> \
  -CHAIN <UCSC-chain> \
  -REJECT <rejected.vcf> \
  -REFERENCE_SEQUENCE <hg38.fa> \
  -LIFTOVER_MIN_MATCH 0.0 \
  -WARN_ON_MISSING_CONTIG true \
  -VALIDATION_STRINGENCY LENIENT
```

`LIFTOVER_MIN_MATCH=0.0` is retained from the original project workflow.
Picard's rejection VCF is retained and can contribute its rejection reason and
attempted locus to a merged record produced by another tool.

## Transanno

Container version: Transanno 0.4.5, built from its tagged source release.

```text
transanno liftvcf \
  --original-assembly <numeric-hg19.fa> \
  --new-assembly <hg38.fa> \
  --chain <transanno-chain> \
  --vcf <numeric-hg19.vcf.gz> \
  --output <raw.vcf> \
  --fail <rejected.vcf>
```

Transanno uses its own chain, built for the numeric-contig hg19 FASTA. The
UCSC chain used by the other three branches is not substituted. The rule does
not pass optional swap, multimapping, or chain-realignment flags; Transanno's
defaults remain in effect. Its failure VCF and command log are retained.

## Common END-only check

Every raw tool output is passed through the same script after the tool has
finished. This is not a fifth coordinate-liftover method:

- tool `CHROM`, `POS`, `ID`, `REF`, and `ALT` are retained;
- target REF is checked against hg38 but never repaired;
- sequence-resolved `END` is derived from the reported target REF span;
- symbolic insertion `END` equals the reported target POS;
- a symbolic interval END is accepted only when the branch's own chain maps
  both source breakpoints to the target interval already reported by the tool;
- failures go to `end_rejected.vcf`.

The separate raw liftover behavior of a tool therefore cannot silently supply
a source-build endpoint to the merged hg38 VCF.

## Tool outputs

| Tool | Successful VCF | Tool-native failure output | Common END failure |
|---|---|---|---|
| CrossMap | `liftover/crossmap/<sample>.hg38.vcf.gz` | `<sample>.hg38.unmap.vcf` | `<sample>.end_rejected.vcf` |
| bcftools | `liftover/bcftools/<sample>.hg38.vcf.gz` | `<sample>.rejected.vcf` | `<sample>.end_rejected.vcf` |
| Picard | `liftover/picard/<sample>.hg38.vcf.gz` | `<sample>.rejected.vcf` | `<sample>.end_rejected.vcf` |
| Transanno | `liftover/transanno/<sample>.hg38.vcf.gz` | `<sample>.rejected.vcf` plus log | `<sample>.end_rejected.vcf` |

The merged file is a coordinate-agreement view. The per-tool VCFs are the
authoritative place to inspect each tool's complete record representation.

## Inspect source-coordinate collisions

`inspect_variant_collisions.py` is a read-only audit script for source records
that share `CHROM`, `POS`, `SVTYPE`, and `END`. It shows the exact REF/ALT
differences, VCF record and decompressed line numbers, reverse-complement
tests, literal sequence change blocks, `QUAL`, `FILTER`, INFO assembly fields,
and all per-sample FORMAT evidence. It labels GT as phased or unphased, reports
`PS`/haplotype tags when present, and explicitly tests for reciprocal `1|0`
versus `0|1` calls in the same sample.

Inspect every collision in the current source VCF:

```text
python3 liftover/inspect_variant_collisions.py \
  results/liftover/source/HG002_dragen_chr22.hg19.vcf.gz
```

Inspect one position and print every allele plus every change block:

```text
python3 liftover/inspect_variant_collisions.py \
  results/liftover/source/HG002_dragen_chr22.hg19.vcf.gz \
  --position 22:49780461 \
  --show-full-alleles \
  --show-full-evidence-sequences \
  --max-blocks 0
```

The reverse-complement checks test literal sequence identity. The change-block
alignment is intended for manual inspection and is not, by itself, evidence
that two records are or are not the same biological event.
Missing phase fields are reported as unavailable; the script never infers
haplotype assignment from `QUAL` or read-support counts alone.
