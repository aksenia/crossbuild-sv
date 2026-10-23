# Plan: annotating the native vs lifted comparison

## Background

HG002 GRCh37 calls are lifted to hg38 and compared with a DRAGEN run of the
same reads on hg38 ([comparison.md](comparison.md)). The cohort has the same
setup and no truth set, so no GIAB benchmark is used.

WG result (on hetzner, in `crossbuild-sv-dev/runs/wg4`;
[pilot_results.md](pilot_results.md)): matched calls are 98.8% at the exact
native position. 1,696 of 17,480 native calls are unmatched:

| diff_reason | n |
| --- | --- |
| no_hg19_sequence (mostly hg38 centromere models) | 541 |
| near_same_type (mostly INS at the same base, 10-30% length difference) | 512 |
| none | 464 |
| other_type_other_build (mostly BND) | 134 |
| filtered_other_build | 45 |

Each step annotates hg38 calls, joins on `native_id` / `source_id` in
`comparison/<s>.links.tsv.gz`, and adds strata to `metrics.tsv` (pooled by
`liftover/aggregate.py`). Prototype outside the Snakefile, verify on the dev
set (14, 21, Y), then whole genome.

## Step 1. Genomic context strata

Goal: give every call one genomic context label, so unmatched calls can be
attributed to a context and compared with matched calls.

### Source

GIAB genome stratifications v3.6, GRCh38:
<https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.6/GRCh38@all/>
(md5s: `GRCh38-genome-stratifications-md5s.txt` in the same directory).

19 files, pooled into 10 flags:

| Flag | Files |
| --- | --- |
| `gaps` | `OtherDifficult/GRCh38_gaps_slop15kb.bed.gz` |
| `satellites` | `LowComplexity/GRCh38_satellites_slop5.bed.gz` |
| `false_dup` | `OtherDifficult/GRCh38_false_duplications_correct_copy.bed.gz`, `OtherDifficult/GRCh38_false_duplications_incorrect_copy.bed.gz`, `OtherDifficult/GRCh38_collapsed_duplication_FP_regions.bed.gz` |
| `segdups` | `SegmentalDuplications/GRCh38_segdups.bed.gz` |
| `XY_difficult` | `XY/GRCh38_chrX_PAR.bed.gz`, `XY/GRCh38_chrY_PAR.bed.gz`, `XY/GRCh38_chrX_XTR.bed.gz`, `XY/GRCh38_chrY_XTR.bed.gz`, `XY/GRCh38_chrX_ampliconic.bed.gz`, `XY/GRCh38_chrY_ampliconic.bed.gz` |
| `hypervariable` | `OtherDifficult/GRCh38_MHC.bed.gz`, `OtherDifficult/GRCh38_KIR.bed.gz`, `OtherDifficult/GRCh38_VDJ.bed.gz` |
| `lowmap_l100` | `Mappability/GRCh38_nonunique_l100_m2_e1.bed.gz` |
| `TR_ge10k` | `LowComplexity/GRCh38_AllTandemRepeats_ge10001bp_slop5.bed.gz` |
| `TR_201to10k` | `LowComplexity/GRCh38_AllTandemRepeats_201to10000bp_slop5.bed.gz` |
| `TR_51to200` | `LowComplexity/GRCh38_AllTandemRepeats_51to200bp_slop5.bed.gz` |

`XY_difficult`: X/Y pseudoautosomal regions, X-transposed regions and
ampliconic regions, i.e. the parts of X and Y duplicated between or within
the sex chromosomes. Not the whole `chrX_nonPAR` / `chrY_nonPAR`.

Not used: `notin*` complements, union files (`Union/`, `lowmappabilityall`,
`allOtherDifficultregions`), finer splits (`segdups_gt10kb`,
`SimpleRepeat_*`, `AllTandemRepeats_le50bp`), homopolymers, GC content,
`GenomeSpecific/` (GIAB samples only), `Ancestry/`, `Functional*/`.

### Priority order

Each call gets exactly one label: the first context in this order that it
hits, else `unique`.

```
gaps > satellites > false_dup > segdups > XY_difficult > hypervariable
     > lowmap_l100 > TR_ge10k > TR_201to10k > TR_51to200 > unique
```

Rule: most structural and build-specific first. Assembly gaps, satellite
arrays and known GRCh38 errors are where GRCh37 and hg38 differ as
sequence; duplications, X/Y and hypervariable loci next; mappability and
tandem repeats last, longest arrays first. A call in a satellite that is
also a tandem repeat is labelled `satellites`. Keep this order in code,
docs and any change to the flag list.

The labels partition the calls, so counts per context add up to the
total. `region=any|none` (CUPs/DISCREPs) stays a separate axis: it
describes liftover instability, not genome context.

### Preprocessing (once per GIAB version)

Script `regions/build_context.sh`. Needs bedtools (>= 2.30), bgzip, tabix.
bedtools is not in the pipeline image: use a bedtools container (e.g.
biocontainers) or add it to the Dockerfile.

General rules:

- BED is 0-based, half-open. Keep it that way; do not add or subtract 1.
- Sort the same way everywhere: `LC_ALL=C sort -k1,1 -k2,2n` (or
  `bedtools sort`). bedtools merge/subtract on unsorted input give wrong
  results without an error.
- Keep only `chr1`-`chr22`, `chrX`, `chrY` (drops alts, decoys, `chrM`).

1. Download and verify:

   ```
   BASE=https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.6/GRCh38@all
   wget -x -nH --cut-dirs=6 "$BASE/LowComplexity/GRCh38_satellites_slop5.bed.gz"   # same for all 19
   wget "$BASE/GRCh38-genome-stratifications-md5s.txt"
   ```

   `--cut-dirs=6` gives `LowComplexity/GRCh38_...`, the same paths as in
   the md5 file. With the 19 paths in `files.txt`:
   `grep -F -f files.txt GRCh38-genome-stratifications-md5s.txt | md5sum -c`
   (19 lines, all `OK`).

2. Pool per flag (one BED per flag, files as in the table above):

   ```
   zcat <files of the flag> | grep -v '^#' | cut -f1-3 \
     | awk '$1 ~ /^chr([0-9]+|X|Y)$/' \
     | LC_ALL=C sort -k1,1 -k2,2n | bedtools merge -i - > flags/<flag>.bed
   ```

   `bedtools merge` joins overlapping and book-ended intervals by default;
   no `-d`.

3. Apply the priority order (each flag keeps only bases not claimed by a
   higher one):

   ```
   : > claimed.bed; : > context.bed
   for f in gaps satellites false_dup segdups XY_difficult hypervariable \
            lowmap_l100 TR_ge10k TR_201to10k TR_51to200; do
     if [ -s claimed.bed ]; then bedtools subtract -a flags/$f.bed -b claimed.bed
     else cat flags/$f.bed; fi | awk -v OFS='\t' -v l=$f '{print $1, $2, $3, l}' >> context.bed
     cat claimed.bed flags/$f.bed | LC_ALL=C sort -k1,1 -k2,2n \
       | bedtools merge -i - > claimed.tmp && mv claimed.tmp claimed.bed
   done
   LC_ALL=C sort -k1,1 -k2,2n context.bed | bgzip > giab_v3.6_context.hg38.bed.gz
   tabix -p bed giab_v3.6_context.hg38.bed.gz
   ```

   Use plain `bedtools subtract` (removes only the overlapping bases), not
   `-A` (drops the whole interval on any overlap).

4. Sizes: bp per label from the partition
   (`awk '{bp[$4] += $3 - $2}'`); genome size = sum of chr1-22, X, Y
   lengths in `hg38.fa.fai`; `unique` = genome size minus all labels.
   Write `giab_v3.6_context_sizes.tsv` (`context`, `bp`, `fraction`).

Checks:

- No overlaps:
  `zcat giab_v3.6_context.hg38.bed.gz | awk '$1 == c && $2 < e {n++} {c = $1; e = $3} END {print n + 0}'`
  prints 0.
- Partition bp equal the bp of the merged union of all 10 flag BEDs
  (`claimed.bed` after the loop).
- Every label present; no contig outside chr1-22, X, Y.

Commit the script, md5 list and sizes table. The BED (several MB) goes next
to the references on hetzner; path in `config.yaml` (`regions:
giab_context`). Keep the 10 pooled flag BEDs (`flags/`) for ad hoc 0/1
columns.

### Script

The steps above are best collected in one script, so the context BED can be
rebuilt the same way for a new GIAB version or by someone else. Suggested
shape:

- `regions/giab_context_flags.tsv` with columns `priority`, `flag`, `file`
  (one row per GIAB file). Keeping the flag list and priority order in one
  small table means the script loops over it instead of hard-coding names,
  and changing the order or adding a flag is a one-line edit.
- `regions/build_context.sh <giab_version> <hg38.fa.fai> <outdir>`, writing
  the context BED (+ `.tbi`), the sizes table and the per-flag BEDs to
  `<outdir>`.
- `set -euo pipefail` at the top, so a failed download or command stops the
  script instead of producing a partial BED.
- Run the checks at the end and exit with an error if one fails; then a
  broken BED never reaches the pipeline.
- Print per flag the number of intervals and bp before and after the
  subtraction. This shows how much each flag loses to higher-priority ones
  and is a quick sanity check of the order.

### Lookup

In `compare_native.py`, by breakpoint, not span:

- DEL/DUP/INV: POS and END; INS: POS.
- Imprecise calls: each breakpoint widened by its CIPOS / CIEND.
- Label: highest-priority context hit by any breakpoint.

Output: `context` column in `links.tsv.gz` and `context=<label>` strata in
`metrics.tsv`. Native calls use their own coordinates; lifted calls their
hg38 placement.

### Report

All shares are counts of calls (not bp), per context label:

- **share of matched calls**: matched calls with this label / all matched
  calls;
- **share of unmatched calls**: native-only (or lifted-only) calls with this
  label / all native-only (or lifted-only) calls;
- **genome fraction**: bp with this label / genome size, from
  `giab_v3.6_context_sizes.tsv` (preprocessing step 4).

Two ratios:

- **unmatched vs matched** = share of unmatched / share of matched. Main
  measure: does this context hold more of the discordant calls than of the
  concordant ones? Matched calls are the right baseline because SVs are not
  spread evenly over the genome (most sit in repeats anyway).
- **unmatched vs genome** = share of unmatched / genome fraction. Secondary:
  how concentrated discordance is relative to the size of the context.

Example (made-up numbers): if 400 of 1,696 native-only calls (23.6%) are
`satellites`, but only 150 of 15,778 matched calls (0.95%), the ratio is
about 25: satellites hold far more discordant than concordant calls. If
satellites cover 3% of the genome, unmatched vs genome is about 8.

Check on the dev set that the counts per context add up to the totals.

## Step 2. REF-changed loci: one allele-state tag

Goal: label, the same way for all tools, whether the hg38 reference carries
the variant allele.

### Background

171 source calls lie where the reference allele changed between builds
(`refchanged` set). 29 match a native call (precision 0.17). Two cases:

- **swap**: the hg19 ALT is the hg38 REF. A 1/1 call is hom-ref on hg38, so
  no native call is expected; a het call becomes the opposite event.
- **third allele**: hg38 has a different allele; a native call is expected.

Tool output at these loci:

| | CrossMap | bcftools | Picard | transanno |
| --- | --- | --- | --- | --- |
| flag | none | `SWAP=1` / `SWAP=-1` | none | `REF_CHANGED` |
| swap, hom (GT set to 0/0) | | 38 | | 44 |
| swap, het | | 16 | | 14 |
| third allele | | 100 | | 81 |
| no flag | 154 (GT kept) | 3 | 46 (GT kept) | 3 |
| rejected | 17 | 14 | 125 | 29 |

bcftools and transanno agree on 115 IDs and conflict on 7. Currently
`compare_native.py` uses the source GT, so lifted 0/0 records count as
lifted-only, and swapped records keep the old SVTYPE.

### Plan

1. One tag per source ID on every tool's record, e.g.
   `INFO/HG38_REF=swap|third_allele`. Rule: A (any tool says swap, 66 IDs),
   B (swap unless a tool says third allele, 59 IDs) or C (check against the
   hg38 FASTA).
2. `compare_native.py`: lifted-only + `HG38_REF=swap` + hom source GT gives
   `diff_reason=hg38_ref_matches_sample`, counted separately.
3. Optional: set SVTYPE of swapped records to match their hg38 alleles.
4. Report the change in `refchanged` metrics.

Data (in `crossbuild-sv-dev/runs/wg4`): `liftover/<tool>/HG002.hg38.raw.vcf.gz`,
`liftover/merged/HG002.refchanged_ids.tsv`, links rows with `set=refchanged`.

## Step 3. Same-site insertions with different length

Goal: test whether the ~230 same-base INS pairs with 10-30% length
difference are one repeat expansion assembled to different lengths.

- Tandem repeat catalog (adotto or UCSC simpleRepeat): motif and copy number.
- Is the length difference a whole number of motifs? Compare inserted
  sequences (Truvari sequence similarity).
- Decide whether "same site, different allele" becomes a `diff_reason`.

## Step 4. AnnotSV

Goal: functional and population context of unmatched calls.

- AnnotSV on native and lifted hg38 calls (genes, exons, gnomAD-SV/DGV
  frequency, known pathogenic/benign SVs, ACMG-like class), joined by ID on
  full rows.
- Do unmatched calls hit exons or disease genes? Are they mostly common? Do
  matched pairs get the same class?
- Large annotation download; run as a separate container step.

## Step 5. Chain context

Goal: relate remaining discordance to local assembly differences.

- Per call: distance to chain gap or block edge, chain score, inverted or
  secondary chain, overlap with hg38 fix patches or ALT loci.
- Reuse `ChainIndex` in `liftover/finalize_sv.py`.
