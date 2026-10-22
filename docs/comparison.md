# Comparison with native hg38 calls

Runs when `native_vcf` is set in the config (`snake/rules/comparison.smk`,
`liftover/compare_native.py`). The native hg38 calls of the same individual
are the reference: a lifted call is judged by how well it reproduces the
native call, not the other way round.

## Steps

1. **prep_native_vcf** — the native VCF gets the same treatment as the
   source: SVTYPE/PASS filter, same contig subset (in chr naming), ID stamp
   `CHROM_POS_SVTYPE_END` (hg38 coords, `_2`, `_3` on collisions),
   `bcftools norm` + `finalize_sv.py` (REF check, END from the record), hg38
   region flags (CUPs, DISCREPs). Records failing QC go to
   `comparison/native/<s>.qc_rejected.vcf`.
2. **lifted_set** — each tool's normalized hg38 VCF is split into
   `main` (source IDs not REF-changed) and `refchanged` (IDs any tool flagged
   as REF-changed, see `refchanged_ids` in `liftover.smk`).
3. **svdb_link** — `svdb --merge --no_intra` per tool × set, native listed
   first, `--overlap 0.8`, `--data_profile sample` (config `comparison:`).
   SVDB is used **only** to decide which native and lifted calls belong
   together; no SVDB coordinates or statistics are reported.
4. **compare_native** — reads the SVDB clusters, looks every member up in the
   normalized VCFs and writes the two tables below.

`tool=union` is not a separate SVDB run; it is derived from the per-tool
links (one SVDB run over the merged set would let a native call join only one
of a source ID's several placements):

- every distinct matched (native, source, placement) from any tool;
- `lifted_only` placements of source IDs that no tool matched;
- `native_only` calls that no tool matched.

Both files carry `sample_id` in every row, so a cohort is their
concatenation (e.g. DuckDB
`read_csv('*/comparison/*.metrics.tsv', delim='\t')`).

## `comparison/<sample>.links.tsv.gz`

One row per (tool, set, link). `.` = not applicable / missing.

| Column | Meaning |
|---|---|
| `sample_id`, `tool`, `set` | `tool` ∈ configured tools + `union`; `set` ∈ `main`, `refchanged` |
| `match_class` | `matched` (native + lifted in one SVDB cluster), `native_only`, `lifted_only`. The refchanged set has no `native_only` rows (all native calls are linked against it, so "unmatched" is not meaningful there) |
| `native_id`, `source_id` | stamped IDs: native hg38, source hg19 |
| `n_links` | matched rows: number of distinct source IDs linked to this native call within (tool, set). >1 = one native call matched several lifted calls |
| `nat_chrom/pos/end/svtype` | native record (normalized) |
| `nat_svlen` | native event size in bp, positive, from the record (allele lengths; SVLEN for symbolic INS; END−POS for symbolic intervals) |
| `nat_gt` | native GT, first sample |
| `nat_precise` | 1 unless INFO/IMPRECISE |
| `nat_cipos`, `nat_ciend` | native CIPOS/CIEND as `lo,hi` (`0,0` if absent) |
| `nat_CUPs`, `nat_DISCREPs` | hg38 region flags of the native call (0/1) |
| `src_chrom/pos/end/svtype/svlen/gt` | hg19 source record |
| `src_CUPs/DISCREPs/CENTEL/SEGDUP` | hg19 region flags of the source call |
| `lift_chrom/pos/end` | the hg38 placement (normalized) of this source ID by this tool |
| `lift_svlen` | event size of the lifted record |
| `tools_same_placement` | tools that put the source ID at exactly this CHROM/POS/END |
| `n_tools_lifted` | tools that lifted this source ID at all (any placement) |
| `n_placements` | distinct CHROM/POS/END placements of this source ID across tools |
| `tools_rejected` | `tool:reason` for tools that did not lift it; reason is the tool's own (CrossMap `Fail(...)`, bcftools/Picard FILTER, transanno `FAILED_REASON`) or `qc:<reason>` from `finalize_sv.py` |
| `lift_CUPs`, `lift_DISCREPs` | hg38 region flags at the lifted placement |
| `exact` | 1 if same chrom and `d_pos == d_end == 0` |
| `d_pos`, `d_end` | lifted − native POS / END (same chrom only) |
| `d_svlen` | `lift_svlen − nat_svlen` |
| `recip_overlap` | overlap / longer interval, DEL/DUP/INV only (`.` for INS) |
| `within_ci` | `d_pos` within native CIPOS and `d_end` within CIEND (for precise calls this equals `exact`) |
| `chrom_match`, `svtype_match` | lifted chrom = native chrom; source SVTYPE = native SVTYPE |
| `gt_match` | unphased, allele-sorted source GT = native GT |
| `closest_tool` | tool(s) whose placement of this source ID is closest to the native call (`|d_pos|+|d_end|`, all tools, ties comma-joined) |
| `svtype`, `size_bin`, `precision` | strata, taken from the native call when present (it is the reference), else from the source call. Size bins: `<50`, `50-100`, `100-500`, `500-1k`, `1k-10k`, `>=10k` |
| `in_region` | 1 if any native, source or lifted region flag is set |

## `comparison/<sample>.metrics.tsv`

Long format: `sample_id, tool, set, stratum, metric, value`. Counts are
reported next to every rate so a cohort can pool numerators and denominators
instead of averaging rates.

### `set=all` — liftover yield (no native involved)

Strata: `all`, `svtype=<T>`. Counted over source IDs.

| tool | metric | meaning |
|---|---|---|
| each tool | `source_ids` | filtered source records |
| | `lifted_ids` | in the tool's normalized hg38 output |
| | `tool_rejected_ids` | in the tool's own reject file (`.unmap` / `rejected.vcf`) |
| | `qc_rejected_ids` | lifted by the tool but failed `finalize_sv.py` QC |
| | `unaccounted_ids` | neither lifted nor rejected — should be 0 |
| | `refchanged_flagged_ids` | IDs this tool flagged as REF-changed |
| | `yield` | `lifted_ids / source_ids` |
| `union` | `source_ids`, `lifted_ids`, `yield` | as above, lifted by any tool |
| | `all_tools_agree_ids` | one placement shared by all tools |
| | `multi_placement_ids` | tools disagree on CHROM/POS/END |
| | `refchanged_ids` | IDs routed to the refchanged set |
| | `all_tools_agree_rate` | `all_tools_agree_ids / lifted_ids` |

### `set=main` / `set=refchanged` — agreement with native

Strata: `all`, `svtype=`, `size=`, `precision=`, `region=any|none`.

| metric | meaning |
|---|---|
| `lifted_ids` | source IDs in the set (matched or lifted_only) |
| `matched_ids` | source IDs linked to a native call |
| `lifted_only_ids` | source IDs with no native link |
| `native_ids`, `native_matched`, `native_only` | native calls; linked; not linked (main only) |
| `precision` | `matched_ids / lifted_ids` |
| `recall` | `native_matched / native_ids` (main only) |
| `exact_ids`, `exact_rate` | matched IDs with `exact=1`; / `matched_ids` |
| `within_ci_ids`, `within_ci_rate` | matched IDs with `within_ci=1`; / `matched_ids` |
| `gt_match_ids`, `gt_match_rate` | matched IDs with `gt_match=1`; / `matched_ids` |
| `median_abs_d_pos`, `median_abs_d_end` | over matched rows on the same chrom |

Caveats: an ID counts as exact/within-CI if any of its matched rows is; a
stratum taken from the native call can differ from the source call's
(e.g. size bin), so `lifted_ids` per stratum is not a partition of the source.
