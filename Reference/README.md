# Local reference assemblies

Reference sequences are required to execute the workflow but are intentionally
not stored in Git. Configure or mount these files at the paths declared in
`snake/config.yaml`:

```text
Reference/hg19/hg19_nochr_fixed.fa
Reference/hg19/hg19_nochr_fixed.fa.fai
Reference/hg19/hg19.fa
Reference/hg19/hg19.fa.fai
Reference/hg38/hg38.fa
Reference/hg38/hg38.fa.fai
Reference/hg38/hg38.dict
```

The hg19 numeric-contig and chr-prefixed FASTAs must describe the same source
assembly. The hg38 FASTA must use chr-prefixed contigs. Reference provenance
and checksums should be recorded when the production input set is finalized.

See `docs/pipeline.md` for the mount layout and execution commands.
