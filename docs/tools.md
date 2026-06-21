# Liftover tools — SV-specific notes

All four tools are included in the `crossbuild-sv` container and run on the
same filtered, ID-stamped hg19 VCF. Their outputs are merged on the
`CHROM/POS/ID` key (ID = `CHROM_POS_SVTYPE_END` from hg19).

---

## CrossMap

**Version tested:** 0.7.x  
**Input:** chr-prefixed hg19 VCF (after `normalize_contigs`)

CrossMap has two subcommands relevant to SVs:

- **`CrossMap vcf`** — what this pipeline uses. Lifts VCF records by mapping
  POS through the UCSC chain file and updates INFO/END (since v0.4.3 for
  sequence-resolved SVs; since v0.7.0 for symbolic alleles like `<DEL>`).
- **`CrossMap region`** — takes BED-format input (not VCF). Treats each span
  as a unified region, calculates the fraction of bases that successfully
  remap, and only outputs regions where this fraction exceeds a threshold
  (default 0.85, controlled by `-r`). Useful as a pre-filter to identify SVs
  that span chain-break regions with low mapping confidence — but requires
  converting VCF coordinates to BED and back.

**SV behaviour (`CrossMap vcf`):**

- Symbolic alleles (`<DEL>`, `<INS>` …): supported since v0.7.0; passed through
  without REF sequence verification — liftover proceeds as long as POS maps.
- Sequence-resolved SVs: REF allele at the new position is verified against
  the hg38 reference by default. Records where REF mismatches go to the
  `.unmap` file. Use `--no-comp-alleles` to skip this check and allow records
  through even when the destination REF differs (useful for DRAGEN-style SVs
  where the padded REF base may not match perfectly).
- INFO/END: remapped via the chain independently from POS.

**Key flag:**

- `--no-comp-alleles`: skip destination REF allele verification. Increases
  liftover yield for sequence-resolved SVs at the cost of not detecting
  allele discordances.

**Known limitations:**

- Records that fail to lift appear in `<sample>.hg38.unmap` (not in the main
  rejected VCF — check this file separately).
- `CrossMap region` is not used in this pipeline but is worth running
  separately to flag SVs with low chain-mapping ratios.

**Pipeline flags used:** none beyond defaults (no `--no-comp-alleles`).

---

## bcftools +liftover

**Version tested:** 1.22 with freeseek liftover plugin  
**Input:** chr-prefixed hg19 VCF

Uses the UCSC chain file and verifies REF alleles against both the source and
destination FASTA. One of few tools confirmed to update INFO/END correctly
(Table 2 of the bcftools/liftover paper, Genovese et al. 2024).

**SV behaviour:**

- Handles both symbolic and sequence-resolved alleles.
- INFO/END: by default, END is **recomputed** from POS + the lifted allele
  length. Alternatively, `--lift-end` maps the END coordinate directly through
  the chain (useful when you want chain-based rather than length-based END
  updating, e.g. for large SVs where allele length is uncertain).
- Rejected records go to `<sample>.rejected.vcf` with a per-record rejection
  reason tag.

**Pipeline flags used:**

- `--no-left-align`: **critical for SVs**. Without this, bcftools re-anchors
  the variant to the leftmost unambiguous position after lifting. For a large
  DEL or INS this shifts POS, which then corrupts INFO/END (END is recomputed
  from the shifted POS). Disabling left-alignment preserves the original
  breakpoint coordinates as faithfully as the chain allows.

**Note on `--lift-end`:** not used in this pipeline; END is recomputed from
the lifted POS. If testing shows END errors for large SVs, switching to
`--lift-end` should be the first thing to try.

---

## Picard LiftoverVcf

**Version tested:** 3.4.0  
**Input:** chr-prefixed hg19 VCF (header patched to VCFv4.3)

Picard performs liftover by walking the chain file and mapping each coordinate.
It is conservative: it rejects records that do not meet its confidence
thresholds rather than producing uncertain liftovers.

**SV behaviour:**

- Symbolic alleles (`<DEL>`, `<INS>` …): carried through without REF sequence
  verification. Picard does not attempt to validate that the reference matches
  at the destination — the symbolic allele is simply propagated to hg38.
- Sequence-resolved SVs: REF is verified at the new position; mismatches cause
  rejection.
- INFO/END: lifted independently from POS by looking up the END coordinate in
  the chain. For SVs that fall entirely within a single chain interval this is
  accurate. For SVs that span a chain interval boundary, POS and END may be
  resolved from different intervals, producing a plausible but lower-confidence
  result.

**Pipeline flags used:**

- `LIFTOVER_MIN_MATCH=0.0`: Picard's default is **1.0** (100% of the variant's
  bases must map to the same chain interval). For large SVs spanning a chain
  break this causes immediate rejection even when most of the variant maps
  cleanly. Setting to 0.0 disables the threshold entirely and lets Picard
  attempt the liftover regardless of interval coverage.
- `VALIDATION_STRINGENCY LENIENT`: prevents hard failures on minor VCF
  formatting quirks.
- Header is patched from VCFv4.4 → VCFv4.3 before Picard sees the file;
  Picard rejects anything declared above v4.3.

**Known limitations:**

- Rejected records go to `<sample>.rejected.vcf` and are annotated with
  `PICARD_REJECT_REASON` in the merged output.

---

## Transanno

**Version tested:** 0.4.5  
**Input:** bare-contig hg19 VCF (numeric, no chr prefix)

Transanno uses a minimap2-generated chain (separate from the UCSC chain used by
the other three tools). It is primarily designed for SNV/indel liftover.

**SV behaviour:**

- Symbolic alleles (`<DEL>`, `<INV>` …): not explicitly documented as
  supported. Transanno is documented to handle asterisk (`*`) and dot (`.`) in
  the ALT column; behaviour with VCF symbolic alleles is unclear and rejection
  rates may be high.
- Sequence-resolved SVs: transanno verifies the REF sequence at the new
  position. Large SVs (e.g. a 500 bp DEL with full sequence in REF) are more
  likely to be rejected due to strict REF verification.
- INFO/END: updated when the record lifts successfully.
- If REF and ALT are swapped (reference flipped at the new position), transanno
  rewrites REF, ALT, AF, and GT accordingly. Use `--noswap` to disable this
  behaviour.
- `MULTIMAP` flag: set when the variant maps to more than one location in hg38.
  Common in segmental duplication regions. Use `--allow-multi-map` to keep
  these records instead of sending them to the fail file.

**Pipeline flags used:** none beyond defaults.

**Known limitations:**

- Expects a minimap2-generated chain in transanno's own format — not the
  standard UCSC `.over.chain`. See `preprocess/` docs for chain build
  instructions.
- `--no-left-align-chain` (note: different from bcftools's `--no-left-align`)
  prevents the provided chain from losing its 1-to-1 mapping properties during
  processing; may be needed for chains built from assemblies with complex
  rearrangements.
- Check `<sample>.log` and `<sample>.rejected.vcf` for rejection details.
- The transanno chain is built against numeric-contig reference FASTAs; the
  pipeline strips chr prefixes before running transanno.

---

## Tool comparison summary

| | CrossMap | bcftools | Picard | Transanno |
| --- | --- | --- | --- | --- |
| Symbolic alleles | ✓ (v0.7.0+) | ✓ | ✓ (carried through) | unclear |
| Sequence-resolved SVs | ✓ | ✓ | ✓ | ✓ (strict REF check) |
| INFO/END update | ✓ | ✓ (recompute or `--lift-end`) | ✓ (per-interval) | ✓ |
| Rejection log | `.unmap` file | rejected VCF | rejected VCF | rejected VCF + log |
| Chain format | UCSC | UCSC | UCSC | minimap2 (transanno) |
| Key SV flag | `--no-comp-alleles` | `--no-left-align` | `LIFTOVER_MIN_MATCH=0.0` | `--noswap`, `--allow-multi-map` |
| Notable alternative | `CrossMap region` (BED) | `--lift-end` for large SVs | default match=1.0 | `--no-left-align-chain` |
