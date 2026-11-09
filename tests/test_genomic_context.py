"""Code-only context checks using synthetic intervals and the small dev fixture.

No downloads, external liftover tools, or full reference genomes are needed.
Run from the repo root: python3 -m unittest tests.test_genomic_context -v
"""
import csv
import gzip
import sys
import tempfile
import unittest
from collections import Counter
from pathlib import Path
from unittest.mock import patch

from liftover import compare_native as cn
from tests.test_compare_native import DATA, TOOLS, run_compare, write_chain


class GenomicContextTest(unittest.TestCase):
    def setUp(self):
        self.priorities = {"satellites": 1, "segdups": 2}
        self.regions = {"chr1": [(100, 200, "satellites"), (200, 300, "segdups")]}

    def record(self, pos, svtype="INS", **info):
        return {"chrom": "chr1", "pos": pos, "info": {"SVTYPE": svtype, **info}}

    def test_coordinates_and_half_open_boundaries(self):
        self.assertEqual(cn.breakpoint_interval(101), (100, 101))
        self.assertEqual(cn.breakpoint_interval(101, (-10, 20)), (90, 121))
        for pos, expected in [(100, "unique"), (101, "satellites"),
                              (200, "satellites"), (201, "segdups"), (301, "unique")]:
            with self.subTest(pos=pos):
                self.assertEqual(cn.call_context(self.record(pos), self.regions,
                                                 self.priorities), expected)

    def test_priority_ci_and_breakpoints_not_span(self):
        label = lambda record: cn.call_context(record, self.regions, self.priorities)
        self.assertEqual(label(self.record(205, CIPOS="-10,10")), "satellites")
        for svtype in ("DEL", "DUP", "INV"):
            with self.subTest(svtype=svtype):
                self.assertEqual(label(self.record(50, svtype, END="350")), "unique")
                self.assertEqual(label(self.record(150, svtype, END="250")), "satellites")
                self.assertEqual(label(self.record(50, svtype, END="205", CIEND="-10,10")),
                                 "satellites")

    def test_unassessed_is_not_unique(self):
        self.assertEqual(cn.call_context(self.record(150, "DEL"), self.regions,
                                         self.priorities), cn.NA)
        self.assertEqual(cn.call_context(self.record(150, CIPOS="bad"), self.regions,
                                         self.priorities), cn.NA)
        outside = {**self.record(150), "chrom": "chrM"}
        self.assertEqual(cn.call_context(outside, self.regions, self.priorities), cn.NA)

    def test_map_rejects_overlap_and_unknown_labels(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "context.bed"
            for contents in ("chr1\t100\t200\tsatellites\nchr1\t199\t300\tsegdups\n",
                             "chr1\t100\t200\tunknown\n"):
                path.write_text(contents)
                with self.assertRaises(ValueError):
                    cn.load_context_intervals(path, self.priorities)

    def test_counts_deduplicate_and_keep_sides_separate(self):
        matched = {
            "tool": "crossmap", "set": "main", "match_class": "matched",
            "native_id": "N1", "source_id": "S1", "lift_chrom": "chr1",
            "lift_pos": 201, "lift_end": 201,
            "native_context": "satellites", "lifted_context": "segdups",
        }
        native_only = {**matched, "match_class": "native_only", "source_id": cn.NA}
        lifted_only = {**matched, "match_class": "lifted_only", "native_id": cn.NA}
        counts = cn.context_counts([matched, matched, native_only, lifted_only])
        self.assertEqual(counts[("crossmap", "main", "native", "satellites", "matched")], 1)
        self.assertEqual(counts[("crossmap", "main", "lifted", "segdups", "matched")], 1)
        self.assertEqual(sum(counts.values()), 2)
        # A different END is a different placement, even with the same POS.
        counts = cn.context_counts([matched, {**lifted_only, "lift_end": 250}])
        self.assertEqual(sum(counts.values()), 3)

    def test_shares_and_both_enrichment_ratios(self):
        counts = Counter({
            ("x", "main", "native", "satellites", "matched"): 5,
            ("x", "main", "native", "unique", "matched"): 95,
            ("x", "main", "native", "satellites", "unmatched"): 20,
            ("x", "main", "native", "unique", "unmatched"): 80,
        })
        summaries = cn.context_summary(counts, self.priorities,
                                       {"satellites": 0.02, "segdups": 0.03, "unique": 0.95})
        satellite = next(row for row in summaries if row["context"] == "satellites")
        self.assertAlmostEqual(satellite["matched_share"], 0.05)
        self.assertAlmostEqual(satellite["unmatched_share"], 0.20)
        self.assertAlmostEqual(satellite["unmatched_vs_matched"], 4)
        self.assertAlmostEqual(satellite["unmatched_vs_genome"], 10)
        self.assertAlmostEqual(sum(row["matched_share"] for row in summaries), 1)
        self.assertAlmostEqual(sum(row["unmatched_share"] for row in summaries), 1)

    def test_missing_baselines_and_refchanged_native(self):
        counts = Counter({("x", "refchanged", "native", "satellites", "matched"): 2})
        summaries = cn.context_summary(counts, self.priorities)
        for row in summaries:
            self.assertEqual(row["unmatched_count"], cn.NA)
            self.assertEqual(row["unmatched_share"], cn.NA)
            self.assertEqual(row["unmatched_vs_matched"], cn.NA)

    def test_union_context_agreement(self):
        builder = object.__new__(cn.Builder)
        builder.tools = ("a", "b")
        builder.lifted_context = {
            ("a", "S", "chr1", 150, 150): "satellites",
            ("b", "S", "chr1", 150, 150): "satellites",
        }
        self.assertEqual(builder.placement_context("union", "S", "chr1", 150, 150),
                         "satellites")
        builder.lifted_context[("b", "S", "chr1", 150, 150)] = "segdups"
        self.assertEqual(builder.placement_context("union", "S", "chr1", 150, 150), cn.NA)

    def test_report_cli_and_existing_metrics(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            baseline = root / "baseline"
            baseline.mkdir()
            baseline_links, baseline_metrics = run_compare(str(baseline))
            chain = root / "fixture.chain"
            write_chain(str(chain))
            flags = root / "flags.tsv"
            flags.write_text("priority\tflag\tfile\n1\tsatellites\ta\n2\tsegdups\tb\n")
            bed = root / "context.bed.gz"
            with gzip.open(bed, "wt") as handle:
                handle.write("chr14\t0\t300000000\tsatellites\n"
                             "chr21\t0\t300000000\tsegdups\n"
                             "chrY\t0\t300000000\tsatellites\n")
            sizes = root / "sizes.tsv"
            sizes.write_text("context\tbp\tfraction\nsatellites\t600000000\t0.6666666667\n"
                             "segdups\t300000000\t0.3333333333\nunique\t0\t0\n")
            links, metrics, report = root / "links.tsv.gz", root / "metrics.tsv", root / "context.tsv"
            f = lambda name: str(Path(DATA) / name)
            argv = [
                "compare_native.py", "--sample", "S", "--tools", ",".join(TOOLS),
                "--source", f("source.vcf"), "--native", f("native.vcf"),
                "--tool-vcf", *(f"{t}={f(t + '.vcf')}" for t in TOOLS),
                "--rejected", *(f"{t}={f(t + '.rejected.vcf')}" for t in TOOLS),
                "--qc-rejected", *(f"{t}={f(t + '.qc_rejected.vcf')}" for t in TOOLS),
                "--lifted-annotated", f("lifted.merged.vcf"), f("lifted.refchanged.vcf"),
                "--refchanged-ids", f("refchanged_ids.tsv"),
                "--svdb", *(f"{t}:{s}={f(f'svdb.{t}.{s}.vcf')}" for t in TOOLS for s in cn.SETS),
                "--raw-source", f("source.vcf"), "--raw-native", f("native.vcf"),
                "--chain", str(chain), "--links", str(links), "--metrics", str(metrics),
                "--context-bed", str(bed), "--context-flags", str(flags),
                "--context-sizes", str(sizes), "--context-report", str(report),
            ]
            with patch.object(sys, "argv", argv):
                cn.main()
            with report.open() as handle:
                report_rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertTrue(report_rows)
            self.assertTrue(all(row["sample_id"] == "S" for row in report_rows))
            self.assertEqual({row["count_unit"] for row in report_rows},
                             {"native_call", "lifted_placement"})
            self.assertTrue(any(row["genome_fraction"] != cn.NA for row in report_rows))
            # All pre-existing metrics retain their values when context is enabled.
            def metric_values(path):
                with open(path) as handle:
                    return {(r["tool"], r["set"], r["stratum"], r["metric"]): r["value"]
                            for r in csv.DictReader(handle, delimiter="\t")}
            original = metric_values(baseline_metrics)
            labelled = metric_values(metrics)
            self.assertEqual(original, {key: labelled[key] for key in original})
            with gzip.open(baseline_links, "rt") as handle:
                original_rows = list(csv.DictReader(handle, delimiter="\t"))
            with gzip.open(links, "rt") as handle:
                labelled_rows = list(csv.DictReader(handle, delimiter="\t"))
            added = {"native_context", "lifted_context", "context"}
            strip = lambda row: {k: v for k, v in row.items() if k not in added}
            self.assertEqual([strip(r) for r in original_rows], [strip(r) for r in labelled_rows])


if __name__ == "__main__":
    unittest.main()
