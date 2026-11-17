import os
import tempfile
import unittest

from scripts.inspect_refchanged import IndexedFasta, allele_relation, classification_rows


class InspectRefchangedTest(unittest.TestCase):
    def test_allele_relation_handles_reverse_complement_and_third_allele(self):
        source = {"ref": "A", "alt": "C"}
        self.assertEqual(allele_relation(source, "G"), "swap")
        self.assertEqual(allele_relation(source, "T"), "reference_allele")
        self.assertEqual(allele_relation(source, "N"), "third_allele")
        self.assertEqual(allele_relation(source, "GGA"), "third_allele")

    def test_classification_table_is_one_row_per_id_and_checks_fasta(self):
        with tempfile.TemporaryDirectory() as temp:
            fasta_path = os.path.join(temp, "hg38.fa")
            with open(fasta_path, "wb") as handle:
                handle.write(b">chr1\nACGT\nACGT\n")
            with open(fasta_path + ".fai", "w") as handle:
                handle.write("chr1\t8\t6\t4\t5\n")

            source = {"id1": [{"ref": "A", "alt": "C"}]}
            flagged_by = {"id1": {"bcftools"}}
            data = {
                tool: {"raw": {}, "retained": {}, "tool_rejected": {}, "qc_rejected": {}}
                for tool in ("crossmap", "bcftools", "picard", "transanno")
            }
            data["bcftools"]["raw"]["id1"] = [{
                "SWAP": "none", "ref": "C", "chrom": "chr1", "pos": "2",
            }]

            fasta = IndexedFasta(fasta_path)
            try:
                rows = classification_rows(["id1"], flagged_by, source, data, fasta)
            finally:
                fasta.close()

            self.assertEqual(len(rows), 1)
            self.assertEqual(rows[0], {
                "source_id": "id1",
                "tools": "bcftools",
                "approach_A": "no_swap_flag",
                "approach_B": "swap",
                "approach_C": "swap",
            })

    def test_approach_c_reports_reference_mismatch(self):
        with tempfile.TemporaryDirectory() as temp:
            fasta_path = os.path.join(temp, "hg38.fa")
            with open(fasta_path, "wb") as handle:
                handle.write(b">chr1\nACGT\n")
            with open(fasta_path + ".fai", "w") as handle:
                handle.write("chr1\t4\t6\t4\t5\n")

            source = {"id1": [{"ref": "A", "alt": "C"}]}
            flagged_by = {"id1": {"bcftools"}}
            data = {
                tool: {"raw": {}, "retained": {}, "tool_rejected": {}, "qc_rejected": {}}
                for tool in ("crossmap", "bcftools", "picard", "transanno")
            }
            data["bcftools"]["raw"]["id1"] = [{
                "SWAP": "1", "ref": "A", "chrom": "chr1", "pos": "2",
            }]

            fasta = IndexedFasta(fasta_path)
            try:
                rows = classification_rows(["id1"], flagged_by, source, data, fasta)
            finally:
                fasta.close()

            self.assertEqual(rows[0]["approach_A"], "swap")
            self.assertEqual(rows[0]["approach_C"], "reference_mismatch")

    def test_approach_c_preserves_disagreement_as_conflict(self):
        with tempfile.TemporaryDirectory() as temp:
            fasta_path = os.path.join(temp, "hg38.fa")
            with open(fasta_path, "wb") as handle:
                handle.write(b">chr1\nACGT\n")
            with open(fasta_path + ".fai", "w") as handle:
                handle.write("chr1\t4\t6\t4\t5\n")

            source = {"id1": [{"ref": "A", "alt": "C"}]}
            flagged_by = {"id1": {"bcftools", "picard"}}
            data = {
                tool: {"raw": {}, "retained": {}, "tool_rejected": {}, "qc_rejected": {}}
                for tool in ("crossmap", "bcftools", "picard", "transanno")
            }
            data["bcftools"]["raw"]["id1"] = [{
                "SWAP": "none", "ref": "C", "chrom": "chr1", "pos": "2",
            }]
            data["picard"]["raw"]["id1"] = [{
                "SWAP": "none", "ref": "A", "chrom": "chr1", "pos": "2",
            }]

            fasta = IndexedFasta(fasta_path)
            try:
                rows = classification_rows(["id1"], flagged_by, source, data, fasta)
            finally:
                fasta.close()

            self.assertEqual(rows[0]["approach_C"], "conflict")


if __name__ == "__main__":
    unittest.main()
