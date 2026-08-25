# crossbuild-sv

`crossbuild-sv` is the preprocessing component of an hg19-to-hg38 structural
variant (SV) comparison project. It filters one hg19 SV VCF, runs four
independent liftover tools, merges tool results that agree on the lifted SV
interval, and adds the required genomic-region flags.

The executable workflow ends with two annotated VCFs. Native-hg38 comparison,
concordance metrics, prioritisation, and cohort aggregation are downstream
tasks and are not implemented in this repository.

## Preprocessing workflow

```text
hg19 SV VCF
    |
    +-- filter SVTYPE and FILTER
    +-- stamp a unique source-derived ID
    +-- annotate source: CUPs, DISCREPs, CENTEL, SEGDUP
    |
    +-- CrossMap --------+
    +-- bcftools --------+-- normalize END only
    +-- Picard ----------+-- tag LIFTOVER_TOOL
    +-- Transanno -------+-- merge matching intervals
                              |
                              +-- annotate hg38: CUPs, DISCREPs
```

All four tools receive the same filtered source records. Their successful
outputs and rejection files remain separate and auditable.

## SV-specific coordinate rule

The liftover tools remain the coordinate-producing methods. Post-processing
does not rewrite their `CHROM`, `POS`, `ID`, `REF`, or `ALT`; it only prevents
an hg19 `END` from being mistaken for an hg38 endpoint:

- sequence-resolved records: `END = POS + length(REF) - 1`;
- symbolic `INS`: `END = POS`;
- symbolic `DEL`, `DUP`, `DUP:TANDEM`, and `INV`: source `POS` and `END` must
  map through one branch-specific chain record to the interval reported by the
  tool. The tool's `CHROM` and `POS` must agree, and only `END` is filled;
- a record that cannot satisfy the applicable rule, or whose reported REF does
  not match hg38, is written to the tool's `end_rejected.vcf`.

No endpoint provenance score, event clustering, comparison threshold, or
consensus call is added.

## Source identity and merge

Each retained source record is assigned:

```text
CHROM_POS_SVTYPE_END_rN
```

The coordinate-based prefix follows the original pipeline. The deterministic
record ordinal is necessary because `CHROM_POS_SVTYPE_END` is not unique in the
pilot input.

Tool records are collapsed only when these fields agree:

```text
CHROM / POS / END / source-derived ID
```

`LIFTOVER_TOOL` contains the tools that produced that interval. Different
target positions or endpoints remain separate. REF/ALT differences are not
used to redefine coordinate agreement; the complete tool-specific
representations remain in the per-tool VCFs.

## Supported SV types

The configured set may contain `DEL`, `DUP`, `DUP:TANDEM`, `INS`, and `INV`.
`BND` and any other unconfigured type are removed at the filter step.

## Repository layout

```text
crossbuild-sv/
├── Dockerfile
├── liftover/
│   ├── merge_tools.py
│   └── normalize_sv_coordinates.py
├── regions/
├── snake/
│   ├── Snakefile
│   ├── config.yaml
│   └── rules/
│       ├── liftover.smk
│       └── annotations.smk
└── docs/
    ├── pipeline.md
    └── tools.md
```

`Reference/` contains local reference files for the current run and is not
copied into the image.

## Build and run

From the repository root:

```bash
docker build -t crossbuild-sv:latest .
```

Edit `snake/config.yaml`, then run with the configured reference, data, and
results mounts. For the current HG002 chr22 pilot:

```bash
docker run --rm \
  -v "$PWD/Reference:/ref:ro" \
  -v "$PWD/results:/results" \
  -v "/Users/shaniaim/Documents/OUS/variant-benchmarking/02_chr22/hg19:/data:ro" \
  crossbuild-sv:latest \
  snakemake \
    --snakefile /app/snake/Snakefile \
    --configfile /app/snake/config.yaml \
    --cores 1 \
    --printshellcmds
```

Use `--dry-run` before execution. The first acceptance run uses one core
because the Picard rule reserves a 6 GB Java heap.

See [docs/pipeline.md](docs/pipeline.md) for the complete file-by-file workflow
and [docs/tools.md](docs/tools.md) for the exact four-tool commands and
limitations.

## Final outputs

| File | Preprocessing result |
|---|---|
| `annotated/source/<sample>.hg19.annotated.vcf.gz` | Filtered, ID-stamped hg19 records with CUPs, DISCREPs, CENTEL, and SEGDUP flags |
| `annotated/merged/<sample>.merged.hg38.annotated.vcf.gz` | Merged hg38 liftover intervals with `LIFTOVER_TOOL`, CUPs, and DISCREPs |

Intermediate per-tool VCFs and rejection VCFs are retained under
`liftover/<tool>/`.

These two annotated VCFs are the preprocessing interface to a later
per-sample native-hg38 comparison. Running that comparison across samples and
aggregating its metrics would form cohort mode, but neither stage belongs to
this preprocessing task.
