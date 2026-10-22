"""merge_tools.py: collapse per CHROM/POS/END/ID, REF-changed routing."""
import os
import subprocess
import sys
import tempfile
import unittest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, 'tests', 'data', 'compare')
HEADER = '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'


def body(path):
    with open(path) as fh:
        return [l for l in fh if not l.startswith('#')]


def merge(lines, *args):
    """Run merge_tools on sorted records; return (main records, refchanged-out records)."""
    def key(l):
        p = l.split('\t')
        return p[0], int(p[1])
    with tempfile.TemporaryDirectory() as d:
        rc = os.path.join(d, 'rc.vcf')
        extra = ['--refchanged-ids', args[0], '--refchanged-out', rc] if args else []
        res = subprocess.run([sys.executable, os.path.join(ROOT, 'liftover', 'merge_tools.py'), *extra],
                             input=HEADER + ''.join(sorted(lines, key=key)),
                             text=True, check=True, capture_output=True)
        rc_recs = body(rc) if args else []
    main = [l.split('\t') for l in res.stdout.splitlines() if not l.startswith('#')]
    return main, [l.split('\t') for l in rc_recs]


def tools(rec):
    return dict(f.split('=', 1) for f in rec[7].split(';') if '=' in f)['LIFTOVER_TOOL']


class MergeTest(unittest.TestCase):
    def test_real_records(self):
        lines = body(os.path.join(DATA, 'bcftools.vcf')) + body(os.path.join(DATA, 'picard.vcf'))
        main, rc = merge(lines, os.path.join(DATA, 'refchanged_ids.tsv'))
        by_id = {}
        for p in main:
            by_id.setdefault(p[2], []).append(p)
        # same placement from both tools -> one record
        self.assertEqual([tools(p) for p in by_id['14_19198898_INS_19198898']], ['bcftools,picard'])
        # Picard rejected the symbolic DUP
        self.assertEqual([tools(p) for p in by_id['Y_22298258_DUP_22303389']], ['bcftools'])
        # REF-changed ID: every tool's record leaves the main set, two placements stay apart
        rc_id = '14_22629738_INS_22629738'
        self.assertNotIn(rc_id, by_id)
        self.assertEqual(sorted(tools(p) for p in rc), ['bcftools', 'picard'])
        self.assertTrue(all(p[2] == rc_id and 'REFCHANGED_BY=bcftools' in p[7] for p in rc))

    def test_interleaved_ids_at_one_pos(self):
        rec = 'chr1\t100\t{}\tA\t<DEL>\t.\tPASS\tEND={};LIFTOVER_TOOL={}\n'
        lines = [rec.format('X', 200, 'transanno'), rec.format('Y', 300, 'crossmap'),
                 rec.format('X', 200, 'crossmap')]
        main, _ = merge(lines)
        self.assertEqual([(p[2], tools(p)) for p in main], [('X', 'crossmap,transanno'), ('Y', 'crossmap')])


if __name__ == '__main__':
    unittest.main()
