import gzip
import tempfile
import unittest
from pathlib import Path

from liftover.inspect_variant_collisions import (
    collision_groups,
    direct_hamming,
    exact_reverse_complement,
    load_records,
    minimal_changed_alt,
    reciprocal_phased_genotypes,
    reverse_complement,
)


VCF_TEXT = """##fileformat=VCFv4.3
##INFO=<ID=END,Number=1,Type=Integer,Description="End">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG002
22\t100\tv1\tA\tAAT\t99\tPASS\tEND=100;SVTYPE=INS;CIGAR=1M2I\tGT:GQ:PS:PR:SR\t1|0:80:10:8,4:3,9
22\t100\tv2\tA\tACG\t98\tPASS\tEND=100;SVTYPE=INS;CIGAR=1M2I\tGT:GQ:PS:PR:SR\t0|1:70:10:6,5:2,8
22\t200\tv3\tT\tTAA\t50\tPASS\tSVTYPE=INS\tGT:GQ:PS:PR:SR\t0/1:40:.:4,2:1,5
"""


class InspectVariantCollisionsTest(unittest.TestCase):
    def test_load_and_group_plain_vcf(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "input.vcf"
            path.write_text(VCF_TEXT)
            records = load_records(path)

        groups = collision_groups(records)
        self.assertEqual(len(records), 3)
        self.assertEqual(len(groups), 1)
        self.assertEqual(groups[0][0], ("22", 100, "INS", 100))
        self.assertEqual([record.variant_id for record in groups[0][1]], ["v1", "v2"])
        self.assertEqual(records[2].end, 200)
        self.assertEqual(records[0].qual, "99")
        self.assertEqual(records[0].filters, "PASS")
        self.assertEqual(records[0].samples[0].name, "HG002")
        self.assertEqual(records[0].samples[0].values["GT"], "1|0")
        self.assertEqual(records[0].samples[0].values["PS"], "10")

    def test_load_gzip_vcf(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "input.vcf.gz"
            with gzip.open(path, "wt") as handle:
                handle.write(VCF_TEXT)
            records = load_records(path)

        self.assertEqual(len(records), 3)
        self.assertEqual(records[0].line_number, 5)
        self.assertEqual(records[0].record_number, 1)

    def test_sequence_checks(self):
        self.assertEqual(reverse_complement("ATCG"), "CGAT")
        self.assertTrue(exact_reverse_complement("ATCG", "CGAT"))
        self.assertFalse(exact_reverse_complement("ATCG", "ATCG"))
        self.assertTrue(reciprocal_phased_genotypes("1|0", "0|1"))
        self.assertFalse(reciprocal_phased_genotypes("0/1", "0/1"))
        self.assertIsNone(reciprocal_phased_genotypes(".", "0|1"))
        self.assertEqual(minimal_changed_alt("A", "AAT"), "AT")
        self.assertEqual(direct_hamming("AAT", "ACG"), 2)
        self.assertIsNone(direct_hamming("A", "AA"))


if __name__ == "__main__":
    unittest.main()
