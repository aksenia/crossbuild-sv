# Region annotation BEDs

Copied into the image at `/app/regions/` and referenced from `snake/config.yaml`.
All files are coordinate-sorted, BGZF-compressed, tabix-indexed BEDs
(`CHROM`, 0-based `FROM`, half-open `TO`).

| File | Build | Contigs | INFO flag |
| --- | --- | --- | --- |
| `cups/FASTA_BED.ALL_GRCh37.novel_CUPs.bed.gz` | hg19 | `1`, `2`, … | `CUPs` |
| `cups/FASTA_BED.ALL_GRCh38.novel_CUPs.bed.gz` | hg38 | `chr1`, … | `CUPs` |
| `discreps/DISCREPS_hg19.sorted.bed.gz` | hg19 | `1`, … | `DISCREPs` |
| `discreps/DISCREPS_hg38.sorted.bed.gz` | hg38 | `chr1`, … | `DISCREPs` |
| `centel/centel_hg19.bed.gz` | hg19 | `1`, … | `CENTEL` |
| `segdups/segdups_hg19.bed.gz` | hg19 | `1`, … | `SEGDUP` |

`header.txt` holds the INFO definitions for all four flags.

Chain files and FASTAs are not stored here; they are mounted at runtime
(see `docs/pipeline.md`).
