# crossbuild-sv preprocessing pipeline

This workflow performs only the preprocessing described by the original
`pipeline.md`: filter one hg19 SV call set, lift it independently with four
tools, merge matching lifted intervals, and annotate the source and merged
VCFs with the required region sets.

It does not compare against native hg38 calls and does not calculate
per-sample or cohort metrics.

## Complete workflow

```text
input hg19 SV VCF
       |
       v
filter configured SVTYPEs
optionally require PASS or missing FILTER
stamp source ID: CHROM_POS_SVTYPE_END_rN
       |
       +-----------------------------------------+
       |                                         |
       v                                         v
rename header and record contigs             Transanno
to chr-prefixed hg19                         numeric hg19
       |                                         |
       +-------------+--------------+------------+
       v             v              v            v
    CrossMap      bcftools        Picard       Transanno
                   +liftover
       |             |              |            |
       +-------------+--------------+------------+
                         |
                         v
                normalize END only
                         |
                         v
                 tag LIFTOVER_TOOL
                         |
                         v
        merge CHROM/POS/END/source-derived ID
                         |
             +-----------+-----------+
             v                       v
   annotate source hg19       annotate merged hg38
   CUPs, DISCREPs,            CUPs, DISCREPs
   CENTEL, SEGDUP
```

## 1. Source preparation

`prep_source_vcf` keeps the configured `sv_types` and, when
`pass_only: true`, records whose FILTER is `PASS` or missing. Every retained
record gets this ID:

```text
CHROM_POS_SVTYPE_END_rN
```

`CHROM_POS_SVTYPE_END` is the identity scheme in the original pipeline.
`rN` is a deterministic retained-record ordinal required to distinguish
records that share those four fields. It has no biological interpretation.

The result is sorted, BGZF-compressed, and indexed:

```text
liftover/source/<sample>.hg19.vcf.gz
```

## 2. Contig naming

CrossMap, bcftools, and Picard use the UCSC chain and receive a chr-prefixed
source VCF. The rename map is generated from every `##contig` header entry,
not only contigs observed in records, so the header cannot contain a mixture of
numeric and chr-prefixed names. `MT` maps to `chrM`.

Transanno uses the numeric-contig hg19 reference and its own chain. If an input
is chr-prefixed, all matching header contigs are converted back; `chrM` maps
to `MT`.

## 3. Four independent liftover branches

All branches start from the same filtered, ID-stamped records.

| Branch | Source input | Chain | Required setting |
|---|---|---|---|
| CrossMap | chr-prefixed VCF | UCSC hg19-to-hg38 | `CrossMap vcf --chromid l` |
| bcftools | chr-prefixed VCF | UCSC hg19-to-hg38 | `+liftover --no-left-align` |
| Picard | chr-prefixed VCF with VCFv4.3 header | UCSC hg19-to-hg38 | `LIFTOVER_MIN_MATCH=0.0`, lenient validation |
| Transanno | numeric-contig VCF | Transanno chain | `transanno liftvcf` |

The workflow does not add optional tool modes or pre-filters. Picard receives a
6 GB Java heap because its tested execution otherwise exceeded the available
default.

Each branch writes its own successful VCF and its tool-native failure output.
No failed record is borrowed from another tool.

## 4. Minimal SV END normalization

SVs require an endpoint in addition to the lifted position. Tool output can
otherwise retain an hg19 END or use different END behavior. The common
post-processing step is deliberately limited:

- it does not rewrite tool `CHROM`, `POS`, `ID`, `REF`, or `ALT`;
- it checks that the reported REF matches the hg38 FASTA, without correcting
  the allele;
- it modifies only `INFO/END`.

Rules:

| Source representation | Target END |
|---|---|
| sequence-resolved allele | `POS + length(REF) - 1` |
| symbolic `INS` | `POS` |
| symbolic `DEL`, `DUP`, `DUP:TANDEM`, or `INV` | map source POS and END through one record in that branch's chain; require the mapped target contig and left endpoint to equal the tool's CHROM/POS; use the mapped right endpoint as END |

A symbolic interval without source END is rejected; `SVLEN` is not silently
substituted. A record that fails target REF validation or the applicable END
rule is written to:

```text
liftover/<tool>/<sample>.end_rejected.vcf
```

with `FILTER=END_LIFTOVER_FAILED` and an `END_LIFTOVER_REASON`. Successful
records receive no extra normalization method or confidence annotation.

## 5. Tool tag and merge

Every successful per-tool record receives:

```text
LIFTOVER_TOOL=<tool>
```

The annotation lookup uses `CHROM`, `POS`, and the source-derived ID.

The merged output collapses records only on:

```text
CHROM / POS / END / source-derived ID
```

This is the original position-and-source merge corrected for interval-based
SVs. An END disagreement is a coordinate disagreement and remains on a
separate line. REF and ALT are not extra agreement criteria because the task
is to compare lifted coordinates; each complete tool representation is still
available in its per-tool VCF. When coordinate-matching records have different
payloads, the merged line uses the first representation in the fixed order
CrossMap, bcftools, Picard, Transanno and unions the tool tags.

The merge also carries Transanno's `MULTIMAP` and `REF_CHANGED` fields when
present. If another tool succeeds for a source ID that Picard rejected, the
merged record may include Picard's rejection reason and attempted locus fields.

## 6. Region annotations

These inputs are coordinate-sorted, three-column BED files
(`CHROM`, zero-based `FROM`, half-open `TO`), BGZF-compressed and tabix
indexed. They are annotation resources, not liftover intervals.

| Build/output | BED resources applied |
|---|---|
| source hg19 | CUPs hg19, DISCREPs hg19, CENTEL hg19, SEGDUP hg19 |
| merged hg38 | CUPs hg38, DISCREPs hg38 |

The hg38 output header contains definitions only for the two fields actually
applied there. The hg19 output contains all four.

## Configuration

The current pilot configuration in `snake/config.yaml` is:

```yaml
sample: HG002_dragen_chr22
input_vcf: /vcf/HG002_dragen.sv.hg19.renamed.chr22.svdb.vcf.gz
native_hg38_vcf: /vcf/HG002_dragen.sv.hg38.renamed.chr22.svdb.vcf.gz
results_dir: /results

sv_types: [DEL, INS]
pass_only: true
liftover_tools: [crossmap, bcftools, picard, transanno]
```

The same file supplies:

- numeric and chr-prefixed hg19 FASTAs;
- chr-prefixed hg38 FASTA and its index/dictionary;
- the UCSC and Transanno chains;
- the six build-specific BED resources;
- the two active Python scripts;
- `docker://crossbuild-sv:latest`.

Changing `sv_types` changes the requested filter but does not add another
analysis mode.

## Build and execute

From the `crossbuild-sv` root:

```bash
docker build -t crossbuild-sv:latest .
```

Dry-run the embedded workflow:

```bash
docker run --rm \
  -v "$PWD/Reference:/ref:ro" \
  -v "$PWD/vcf:/vcf:ro" \
  -v "$PWD/results:/results" \
  crossbuild-sv:latest \
  snakemake \
    --snakefile /app/snake/Snakefile \
    --configfile /app/snake/config.yaml \
    --cores 1 \
    --dry-run \
    --printshellcmds
```

Remove `--dry-run` to execute. The output mount must be writable by the
container user.

## Required references

- numeric-contig hg19 FASTA plus `.fai`;
- chr-prefixed hg19 FASTA plus `.fai`;
- chr-prefixed hg38 FASTA plus `.fai` and `.dict`;
- UCSC hg19-to-hg38 chain for CrossMap, bcftools, and Picard;
- Transanno hg19-to-hg38 chain for Transanno;
- sorted, BGZF-compressed, tabix-indexed BED resources listed above.

The UCSC chain and Transanno chain are not interchangeable.

## Outputs

Final preprocessing products:

| File | Contents |
|---|---|
| `annotated/source/<sample>.hg19.annotated.vcf.gz` | filtered and ID-stamped source with four hg19 region flags |
| `annotated/merged/<sample>.merged.hg38.annotated.vcf.gz` | merged interval results with tool tags and two hg38 region flags |

Auditable intermediates:

| File | Contents |
|---|---|
| `liftover/<tool>/<sample>.hg38.vcf.gz` | successful END-normalized output for one tool |
| `liftover/<tool>/<sample>.rejected.vcf` | tool-native failures; CrossMap uses `.hg38.unmap.vcf` |
| `liftover/<tool>/<sample>.end_rejected.vcf` | records excluded only by the explicit common END/REF contract |
| `liftover/transanno/<sample>.log` | Transanno command log |

## Acceptance criteria

Preprocessing is complete when:

1. the filtered source count and SVTYPE/FILTER selection match the config;
2. every retained source ID is unique and is preserved by all tools;
3. all four tool rules execute and retain their independent rejection outputs;
4. normalization changes only END on accepted records;
5. every accepted record has integer `END >= POS` and target-matching REF;
6. merged tool tags exactly reflect agreement on CHROM/POS/END/source ID;
7. the hg19 and hg38 region flags are limited to the sets listed above;
8. both final VCFs are BGZF-compressed and index successfully.

Native-hg38 matching, concordance definitions, thresholds, result
interpretation, and cohort summaries begin only after this boundary.
