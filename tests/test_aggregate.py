"""aggregate.py: pooled counts, rates recomputed from them, n_samples."""
import csv
import os
import subprocess
import sys
import tempfile
import unittest

from liftover.aggregate import RATES
from tests.test_compare_native import run_compare

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
HEADER = 'sample_id\ttool\tset\tstratum\tmetric\tvalue\n'


def aggregate(paths, out):
    subprocess.run([sys.executable, os.path.join(ROOT, 'liftover', 'aggregate.py'),
                    '--output', out, *paths], check=True, capture_output=True)
    with open(out) as fh:
        return list(csv.DictReader(fh, delimiter='\t'))


def read(path):
    with open(path) as fh:
        return {(r['tool'], r['set'], r['stratum'], r['metric']): r['value']
                for r in csv.DictReader(fh, delimiter='\t')}


class AggregateTest(unittest.TestCase):
    def test_two_copies(self):
        with tempfile.TemporaryDirectory() as d:
            paths = []
            for s in ('A', 'B'):
                os.mkdir(os.path.join(d, s))
                paths.append(run_compare(os.path.join(d, s), sample=s)[1])
            one = read(paths[0])
            rows = aggregate(paths, os.path.join(d, 'cohort.tsv'))
        self.assertEqual(len(rows), len(one))
        for r in rows:
            key = (r['tool'], r['set'], r['stratum'], r['metric'])
            self.assertEqual((r['sample_id'], r['n_samples']), ('cohort', '2'))
            v, ref = r['value'], one[key]
            if ref == '.':
                self.assertEqual(v, '.', key)
            elif r['metric'] in RATES or r['metric'].startswith('median_'):
                self.assertEqual(float(v), float(ref), key)
            else:
                self.assertEqual(int(v), 2 * int(ref), key)

    def test_rates_from_pooled_counts(self):
        with tempfile.TemporaryDirectory() as d:
            paths = []
            for s, matched, lifted in (('A', 1, 2), ('B', 3, 4)):
                paths.append(os.path.join(d, s + '.tsv'))
                with open(paths[-1], 'w') as fh:
                    fh.write(HEADER + ''.join(f'{s}\tx\tmain\tall\t{m}\t{v}\n' for m, v in (
                        ('matched_ids', matched), ('lifted_ids', lifted), ('precision', matched / lifted))))
            rows = {r['metric']: r['value'] for r in aggregate(paths, os.path.join(d, 'c.tsv'))}
        self.assertEqual(rows['precision'], '0.6667')   # 4/6, not mean(0.5, 0.75)


if __name__ == '__main__':
    unittest.main()
