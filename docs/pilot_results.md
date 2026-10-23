# Pilot results: HG002 whole genome

Status of the pipeline and numbers from the first whole-genome run with the
native comparison. Metric definitions: [comparison.md](comparison.md).

## Setup

| | |
|---|---|
| Sample | HG002 (GIAB) |
| Source calls | DRAGEN 4.4.4 SV, GRCh37 (`human_g1k_v37_decoy`), 21,659 records |
| Native calls | DRAGEN 4.4.4 SV, same FASTQ, hg38 (`HS_hg38_v5`), 24,221 records |
| Filter | PASS, SVTYPE DEL/DUP/INS/INV (BND out of scope); HG002 has no PASS INV |
| Chains | UCSC hg19ToHg38 (CrossMap, bcftools, Picard); transanno minimap2 chain |
| Tools | CrossMap 0.7.3, bcftools 1.22 +liftover, Picard 3.4.0, transanno 0.4.5 |
| Linking | SVDB 2.12.0 (`--overlap 0.8`, profile `sample`), made 1:1 in compare_native.py |
| Pipeline | main at `f386e54`; run on hetzner, in `crossbuild-sv-dev/runs/wg4`, about 2 min |

After filtering: 16,930 source calls, 17,480 native calls.

**Union** = the four tools combined, derived from the per-tool links (no
separate SVDB run). A source call is lifted if any tool lifted it; a native
call is matched if any tool's placement links to it; native-only and
lifted-only mean no tool matched the call.

## Liftover yield (all source calls)

| Tool | Lifted | Yield | Rejected by tool | Rejected by QC | Unaccounted |
|---|---|---|---|---|---|
| CrossMap | 16,890 | 0.998 | 34 | 6 | 0 |
| bcftools | 16,853 | 0.996 | 71 | 6 | 0 |
| Picard | 16,701 | 0.987 | 229 | 0 | 0 |
| transanno | 16,824 | 0.994 | 101 | 5 | 0 |
| Union | 16,910 | 0.999 | | | |

- All four tools agree on one hg38 placement for 90.4% of lifted calls.
- 171 source calls have a reference allele that changed between builds
  (flagged by at least one tool); they form the separate `refchanged` set.

## Agreement with native calls (main set)

| Tool | Precision | Recall | Exact position | Within native CI | GT match |
|---|---|---|---|---|---|
| CrossMap | 0.941 | 0.902 | 0.913 | 0.925 | 0.998 |
| bcftools | 0.945 | 0.903 | 0.988 | 0.999 | 0.998 |
| Picard | 0.944 | 0.900 | 0.989 | 0.999 | 0.998 |
| transanno | 0.945 | 0.902 | 0.988 | 0.999 | 0.998 |
| Union | 0.943 | 0.903 | 0.988 | 0.999 | 0.998 |

CrossMap's lower exact rate comes from truncating the REF of DRAGEN complex
insertions to one base (see [tools.md](tools.md)).

`refchanged` set (union): 171 lifted, 29 matched to a native call
(precision 0.17), all 29 at the exact position.

### By SV type

| SVTYPE | Native calls | Metric | CrossMap | bcftools | Picard | transanno | Union |
| --- | --- | --- | --- | --- | --- | --- | --- |
| DEL | 8,041 | precision | 0.963 | 0.964 | 0.964 | 0.965 | 0.963 |
| | | recall | 0.924 | 0.924 | 0.924 | 0.924 | 0.924 |
| | | exact position | 0.979 | 0.979 | 0.979 | 0.979 | 0.979 |
| INS | 9,252 | precision | 0.924 | 0.930 | 0.927 | 0.930 | 0.927 |
| | | recall | 0.895 | 0.897 | 0.897 | 0.897 | 0.898 |
| | | exact position | 0.856 | 0.998 | 0.998 | 0.998 | 0.998 |
| DUP | 187 | precision | 0.662 | 0.662 | | 0.682 | 0.662 |
| | | recall | 0.241 | 0.241 | 0.000 | 0.241 | 0.241 |
| | | exact position | 0.622 | 0.622 | | 0.622 | 0.622 |

Picard rejects all symbolic DUPs (`IndelStraddlesMultipleIntervals`).

### Recall by type and size

| SVTYPE | Size | Native calls | CrossMap | bcftools | Picard | transanno | Union |
| --- | --- | --- | --- | --- | --- | --- | --- |
| DEL | <50 | 2,111 | 0.949 | 0.949 | 0.949 | 0.950 | 0.950 |
| | 50-100 | 2,120 | 0.946 | 0.946 | 0.946 | 0.946 | 0.946 |
| | 100-500 | 2,588 | 0.934 | 0.934 | 0.934 | 0.933 | 0.934 |
| | 500-1k | 311 | 0.820 | 0.820 | 0.820 | 0.823 | 0.823 |
| | 1k-10k | 812 | 0.858 | 0.858 | 0.857 | 0.855 | 0.858 |
| | >=10k | 99 | 0.545 | 0.535 | 0.535 | 0.535 | 0.545 |
| INS | <50 | 2,018 | 0.946 | 0.951 | 0.951 | 0.950 | 0.952 |
| | 50-100 | 2,462 | 0.929 | 0.931 | 0.931 | 0.931 | 0.932 |
| | 100-500 | 3,923 | 0.882 | 0.882 | 0.882 | 0.882 | 0.882 |
| | 500-1k | 770 | 0.743 | 0.742 | 0.742 | 0.742 | 0.743 |
| | 1k-10k | 79 | 0.696 | 0.696 | 0.696 | 0.683 | 0.696 |
| DUP | 1k-10k | 169 | 0.237 | 0.237 | 0.000 | 0.237 | 0.237 |
| | >=10k | 18 | 0.278 | 0.278 | 0.000 | 0.278 | 0.278 |

## Position offsets of matched calls

Matched calls by `max(|d_pos|, |d_end|)` in bp:

| SVTYPE | Tool | Matched | 0 | 1-10 | 11-100 | 101-1k | >1k | other chrom | Within CI |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| all | CrossMap | 15,753 | 14,385 | 1,250 | 100 | 18 | 0 | 0 | 0.924 |
| | bcftools | 15,775 | 15,581 | 153 | 25 | 16 | 0 | 0 | 0.999 |
| | Picard | 15,729 | 15,553 | 142 | 22 | 12 | 0 | 0 | 0.999 |
| | transanno | 15,771 | 15,578 | 152 | 25 | 16 | 0 | 0 | 0.999 |
| | Union | 15,778 | 15,582 | 153 | 25 | 18 | 0 | 0 | 0.999 |
| DEL | CrossMap | 7,432 | 7,275 | 126 | 17 | 14 | 0 | 0 | 0.999 |
| | bcftools | 7,431 | 7,275 | 126 | 17 | 13 | 0 | 0 | 0.999 |
| | Picard | 7,430 | 7,275 | 126 | 17 | 12 | 0 | 0 | 0.999 |
| | transanno | 7,429 | 7,274 | 125 | 17 | 13 | 0 | 0 | 0.999 |
| | Union | 7,434 | 7,276 | 126 | 17 | 15 | 0 | 0 | 0.999 |
| INS | CrossMap | 8,276 | 7,082 | 1,113 | 80 | 1 | 0 | 0 | 0.858 |
| | bcftools | 8,299 | 8,278 | 16 | 5 | 0 | 0 | 0 | 0.999 |
| | Picard | 8,299 | 8,278 | 16 | 5 | 0 | 0 | 0 | 0.999 |
| | transanno | 8,297 | 8,276 | 16 | 5 | 0 | 0 | 0 | 0.999 |
| | Union | 8,299 | 8,278 | 16 | 5 | 0 | 0 | 0 | 0.999 |
| DUP | CrossMap | 45 | 28 | 11 | 3 | 3 | 0 | 0 | 0.956 |
| | bcftools | 45 | 28 | 11 | 3 | 3 | 0 | 0 | 0.956 |
| | transanno | 45 | 28 | 11 | 3 | 3 | 0 | 0 | 0.956 |
| | Union | 45 | 28 | 11 | 3 | 3 | 0 | 0 | 0.956 |

- DEL and DUP offsets are the same for every tool: they are differences
  between the two DRAGEN runs (mostly imprecise calls, 120 of 153 bcftools
  1-10 bp offsets are imprecise DELs), not liftover shifts.
- The only tool effect is CrossMap on INS: 1,113 at 1-10 bp and 80 at
  11-100 bp, from truncating the REF of complex insertions.

Offsets of unmatched calls to the nearest same-type call of the other build
(`near_same_type` only, within 1 kb):

| SVTYPE | Tool | Native-only near | 0 | 1-10 | 11-100 | 101-1k | Lifted-only near | 0 | 1-10 | 11-100 | 101-1k |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| all | CrossMap | 535 | 257 | 99 | 69 | 110 | 459 | 277 | 54 | 67 | 61 |
| | bcftools | 519 | 262 | 81 | 66 | 110 | 432 | 251 | 65 | 57 | 59 |
| | Picard | 564 | 293 | 91 | 67 | 113 | 433 | 251 | 66 | 57 | 59 |
| | transanno | 523 | 264 | 82 | 66 | 111 | 430 | 250 | 66 | 56 | 58 |
| | Union | 512 | 257 | 79 | 66 | 110 | 436 | 278 | 40 | 57 | 61 |
| DEL | CrossMap | 103 | 25 | 24 | 10 | 44 | 58 | 13 | 11 | 12 | 22 |
| | bcftools | 104 | 25 | 25 | 10 | 44 | 57 | 13 | 10 | 12 | 22 |
| | Picard | 104 | 25 | 25 | 10 | 44 | 57 | 13 | 10 | 12 | 22 |
| | transanno | 106 | 25 | 26 | 10 | 45 | 57 | 13 | 10 | 12 | 22 |
| | Union | 103 | 25 | 24 | 10 | 44 | 58 | 13 | 11 | 12 | 22 |
| INS | CrossMap | 430 | 232 | 75 | 58 | 65 | 399 | 264 | 43 | 55 | 37 |
| | bcftools | 413 | 237 | 56 | 55 | 65 | 373 | 238 | 55 | 45 | 35 |
| | Picard | 413 | 237 | 56 | 55 | 65 | 376 | 238 | 56 | 45 | 37 |
| | transanno | 415 | 239 | 56 | 55 | 65 | 372 | 237 | 56 | 44 | 35 |
| | Union | 407 | 232 | 55 | 55 | 65 | 376 | 265 | 29 | 45 | 37 |
| DUP | CrossMap | 2 | 0 | 0 | 1 | 1 | 2 | 0 | 0 | 0 | 2 |
| | bcftools | 2 | 0 | 0 | 1 | 1 | 2 | 0 | 0 | 0 | 2 |
| | Picard | 47 | 31 | 10 | 2 | 4 | 0 | 0 | 0 | 0 | 0 |
| | transanno | 2 | 0 | 0 | 1 | 1 | 1 | 0 | 0 | 0 | 1 |
| | Union | 2 | 0 | 0 | 1 | 1 | 2 | 0 | 0 | 0 | 2 |

Picard's 47 native-only DUPs near a same-type hg19 call are the DUPs it
rejected.

## Unmatched calls (main set)

`diff_reason`: what the other build's unfiltered VCF has within 1 kb.

| Native-only | CrossMap | bcftools | Picard | transanno | Union |
| --- | --- | --- | --- | --- | --- |
| near_same_type | 535 | 519 | 564 | 523 | 512 |
| filtered_other_build | 45 | 45 | 45 | 45 | 45 |
| other_type_other_build | 134 | 134 | 134 | 134 | 134 |
| no_hg19_sequence | 542 | 542 | 542 | 541 | 541 |
| none | 465 | 465 | 466 | 466 | 464 |
| **Total** | **1,721** | **1,705** | **1,751** | **1,709** | **1,696** |

| Lifted-only | CrossMap | bcftools | Picard | transanno | Union |
| --- | --- | --- | --- | --- | --- |
| near_same_type | 459 | 432 | 433 | 430 | 436 |
| filtered_other_build | 55 | 53 | 52 | 52 | 55 |
| other_type_other_build | 94 | 91 | 89 | 90 | 95 |
| none | 375 | 345 | 352 | 339 | 375 |
| **Total** | **983** | **921** | **926** | **911** | **961** |

Native-only reasons other than `near_same_type` are identical across tools:
they depend on the hg19 calls, not on liftover. Tools differ only in
`near_same_type`, i.e. calls the tool did not lift or link: Picard +45 (the
rejected DUPs), CrossMap +17 INS compared with bcftools.

Observations (union) from the links table:

- `no_hg19_sequence`: 284 of 541 lie in 29 clusters of 5 or more calls, at
  hg38 centromere and pericentromere models (chr1, 3, 7, 10, 12, 16, 17, 20,
  22, Y). Includes 92 of the 142 native-only DUPs.
- `near_same_type`: 407 of 512 native-only are INS; 232 have an hg19 INS at
  the same base, mostly with a 10-30% different length (SVDB's `sample`
  profile requires SVLEN ratio >= 0.9). Native-only DELs with a same-type
  call 101 bp to 1 kb away mostly differ in size (ratio 0.04-0.8).
- `none`: lower median QUAL than matched calls (588 vs 707), more often
  heterozygous (58% vs 49%) and in flagged regions (21% vs 5%).
- `other_type_other_build`: mostly BND in the other build; 78% heterozygous.
- `filtered_other_build`: median QUAL 322.

## Open points

- No hg38 centromere flag yet; `no_hg19_sequence` is the proxy.
- INS at the same site with a different length are counted as unmatched.
- Next steps: [plan.md](plan.md).
