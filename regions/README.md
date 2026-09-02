# External chain and region resources

The exact chain files and indexed BED payloads consumed by the workflow are
stored in Git and copied into `/app/regions` by the Docker build. Larger raw,
duplicate, and intermediate resources are ignored. The runtime paths and roles
are defined in `snake/config.yaml` and `docs/pipeline.md`.

Expected liftover inputs:

```text
regions/chain/hg19ToHg38.over.chain
regions/chain/transanno.hg19tohg38.chain
```

Expected annotation inputs:

```text
regions/cups/FASTA_BED.ALL_GRCh37.novel_CUPs.bed.gz
regions/cups/FASTA_BED.ALL_GRCh37.novel_CUPs.bed.gz.tbi
regions/cups/FASTA_BED.ALL_GRCh38.novel_CUPs.bed.gz
regions/cups/FASTA_BED.ALL_GRCh38.novel_CUPs.bed.gz.tbi
regions/discreps/DISCREPS_hg19.sorted.bed.gz
regions/discreps/DISCREPS_hg19.sorted.bed.gz.tbi
regions/discreps/DISCREPS_hg38.sorted.bed.gz
regions/discreps/DISCREPS_hg38.sorted.bed.gz.tbi
regions/centel/centel_hg19.bed.gz
regions/centel/centel_hg19.bed.gz.tbi
regions/segdups/segdups_hg19.bed.gz
regions/segdups/segdups_hg19.bed.gz.tbi
```

BED inputs must be coordinate-sorted, BGZF-compressed, and tabix-indexed. The
UCSC chain is used by CrossMap, bcftools, and Picard; Transanno uses its own
chain. The two chain files are not interchangeable.

`regions/header.txt` is also versioned because it supplies the INFO header
definitions used by both annotation rules. Dataset source URLs, versions, and
checksums should be added here when the production input set is finalized.
