"""
compare_native.py on a fixture cut from a dev run (HG002 chr14/21/Y, tools
bcftools + picard). Records: sequence DEL and INS, symbolic DUP (Picard
rejects it), a REF-changed ID, a native-only call, a 1:N SVDB cluster (one
native INS + two near-duplicate source INS). Hand-made: Picard's record of
14_19259092_DEL_19259140 moved to its qc_rejected file (dev data has no QC
rejections). The chain is built here: one 2 kb block around each lifted
record, so hg38 loci away from them have no hg19 sequence.
"""
import csv
import gzip
import os
import subprocess
import sys
import tempfile
import unittest

from liftover import compare_native as cn

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, 'tests', 'data', 'compare')
TOOLS = ('bcftools', 'picard')


def write_chain(path):
    src = {r['id']: r for r in cn.records(os.path.join(DATA, 'source.vcf'))}
    with open(path, 'w') as fh:
        for i, r in enumerate(cn.records(os.path.join(DATA, 'bcftools.vcf'))):
            s, t = src[r['id']]['pos'] - 1001, r['pos'] - 1001
            fh.write(f"chain 1 chr{src[r['id']]['chrom']} 300000000 + {s} {s + 2000} "
                     f"{r['chrom']} 300000000 + {t} {t + 2000} {i}\n2000\n\n")


def run_compare(out_dir, sample='S'):
    """Run compare_native.py on the fixture; return (links, metrics) paths."""
    f = lambda name: os.path.join(DATA, name)
    chain = os.path.join(out_dir, 'fixture.chain')
    write_chain(chain)
    links, metrics = os.path.join(out_dir, 'links.tsv.gz'), os.path.join(out_dir, 'metrics.tsv')
    subprocess.run([
        sys.executable, os.path.join(ROOT, 'liftover', 'compare_native.py'),
        '--sample', sample, '--tools', ','.join(TOOLS),
        '--source', f('source.vcf'), '--native', f('native.vcf'),
        '--tool-vcf', *(f'{t}={f(t + ".vcf")}' for t in TOOLS),
        '--rejected', *(f'{t}={f(t + ".rejected.vcf")}' for t in TOOLS),
        '--qc-rejected', *(f'{t}={f(t + ".qc_rejected.vcf")}' for t in TOOLS),
        '--lifted-annotated', f('lifted.merged.vcf'), f('lifted.refchanged.vcf'),
        '--refchanged-ids', f('refchanged_ids.tsv'),
        '--svdb', *(f'{t}:{s}={f(f"svdb.{t}.{s}.vcf")}' for t in TOOLS for s in cn.SETS),
        '--raw-source', f('source.vcf'), '--raw-native', f('native.vcf'), '--chain', chain,
        '--links', links, '--metrics', metrics], check=True, capture_output=True)
    return links, metrics


class CompareNativeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        with tempfile.TemporaryDirectory() as d:
            links, metrics = run_compare(d)
            with gzip.open(links, 'rt') as fh:
                cls.links = list(csv.DictReader(fh, delimiter='\t'))
            with open(metrics) as fh:
                cls.metrics = {(r['tool'], r['set'], r['stratum'], r['metric']): r['value']
                               for r in csv.DictReader(fh, delimiter='\t')}

    def rows(self, tool, set_='main', **kw):
        return [r for r in self.links if r['tool'] == tool and r['set'] == set_
                and all(r[k] == v for k, v in kw.items())]

    def test_svdb_parsing(self):
        clusters = list(cn.svdb_clusters(os.path.join(DATA, 'svdb.bcftools.main.vcf')))
        multi = [c for c in clusters if len(c['lifted']) > 1]
        self.assertEqual(len(multi), 1)
        self.assertEqual(multi[0]['native'], [('chr21_44918935_INS_44918935', 'chr21', 44918935)])
        self.assertEqual(sorted(i for i, _, _ in multi[0]['lifted']),
                         ['21_46338849_INS_46338850', '21_46338850_INS_46338850'])
        native_only = [c for c in clusters if not c['lifted']]
        self.assertEqual([c['native'][0][0] for c in native_only], ['chr14_20236451_INS_20236451'])

    def test_one_to_one(self):
        m = self.rows('bcftools', native_id='chr21_44918935_INS_44918935')
        self.assertEqual([(r['match_class'], r['source_id'], r['n_links']) for r in m],
                         [('matched', '21_46338849_INS_46338850', '2')])
        demoted = self.rows('bcftools', source_id='21_46338850_INS_46338850')
        self.assertEqual([r['match_class'] for r in demoted], ['lifted_only'])

    def test_rows(self):
        dup = self.rows('bcftools', source_id='Y_22298258_DUP_22303389')[0]
        self.assertEqual((dup['match_class'], dup['exact'], dup['d_pos'], dup['within_ci']),
                         ('matched', '0', '-1', '1'))
        self.assertEqual(dup['tools_rejected'], 'picard:IndelStraddlesMultipleIntervals')
        qc = self.rows('bcftools', source_id='14_19259092_DEL_19259140')[0]
        self.assertEqual((qc['match_class'], qc['tools_rejected']), ('lifted_only', 'picard:qc:REF_MISMATCH'))
        self.assertEqual(self.rows('bcftools', 'refchanged', match_class='native_only'), [])

    def test_union(self):
        # DUP: Picard has it native_only, bcftools matched -> union matched only
        dup = self.rows('union', native_id='chrY_20136373_DUP_20141502')
        self.assertEqual([r['match_class'] for r in dup], ['matched'])
        self.assertEqual(sorted(r['match_class'] for r in self.rows('union')),
                         ['lifted_only'] * 2 + ['matched'] * 4 + ['native_only'])
        # the REF-changed ID has a different END per tool: two placements
        self.assertEqual(len(self.rows('union', 'refchanged')), 2)

    def test_diff_reason(self):
        why = lambda r: (r['diff_reason'], r['nearest_dist'])
        # demoted 1:N partner: the same native INS is right there
        self.assertEqual(why(self.rows('bcftools', source_id='21_46338850_INS_46338850')[0]),
                         ('near_same_type', '0'))
        self.assertEqual(why(self.rows('bcftools', source_id='14_19259092_DEL_19259140')[0]), ('none', '.'))
        # Picard rejected the DUP; the hg19 call lifts 1 bp from the native one
        self.assertEqual(why(self.rows('picard', native_id='chrY_20136373_DUP_20141502')[0]),
                         ('near_same_type', '1'))
        self.assertEqual(why(self.rows('union', native_id='chr14_20236451_INS_20236451')[0]),
                         ('no_hg19_sequence', '.'))
        self.assertTrue(all(r['diff_reason'] == '.' for r in self.links if r['match_class'] == 'matched'))

    def test_other_build(self):
        ob = cn.OtherBuild([('chr1', 1000, 'PASS', 'BND', 0), ('chr1', 1500, 'MinQUAL', 'DEL', 100),
                            ('chr1', 5000, 'PASS', 'DEL', 80)])
        self.assertEqual(ob.explain('chr1', 1200, 'DEL', 100), ('filtered_other_build', '.', '.'))
        self.assertEqual(ob.explain('chr1', 1200, 'INS', 100), ('other_type_other_build', '.', '.'))
        self.assertEqual(ob.explain('chr1', 4500, 'DEL', 100), ('near_same_type', 500, 0.8))
        self.assertEqual(ob.explain('chr2', 1000, 'DEL', 100), ('none', '.', '.'))

    def test_metrics(self):
        m = self.metrics
        expect = {
            ('bcftools', 'main', 'all', 'matched_ids'): '4',
            ('bcftools', 'main', 'all', 'lifted_ids'): '6',
            ('bcftools', 'main', 'all', 'precision'): '0.6667',
            ('bcftools', 'main', 'all', 'recall'): '0.8',
            ('bcftools', 'main', 'all', 'exact_rate'): '0.75',
            ('picard', 'main', 'all', 'recall'): '0.6',
            ('picard', 'all', 'all', 'tool_rejected_ids'): '1',
            ('picard', 'all', 'all', 'qc_rejected_ids'): '1',
            ('picard', 'all', 'all', 'unaccounted_ids'): '0',
            ('picard', 'all', 'all', 'yield'): '0.7143',
            ('union', 'all', 'all', 'all_tools_agree_ids'): '4',
            ('union', 'main', 'all', 'recall'): '0.8',
            ('bcftools', 'main', 'svtype=DEL|size=50-100', 'matched_ids'): '1',
            ('bcftools', 'main', 'all', 'matched_offset_0'): '3',
            ('bcftools', 'main', 'all', 'matched_offset_1-10'): '1',
            ('bcftools', 'main', 'all', 'lifted_only_near_same_type'): '1',
            ('bcftools', 'main', 'all', 'lifted_only_none'): '1',
            ('bcftools', 'main', 'all', 'native_only_no_hg19_sequence'): '1',
            ('picard', 'main', 'all', 'native_only_near_same_type'): '1',
            ('picard', 'main', 'all', 'native_only_near_offset_1-10'): '1',
        }
        self.assertEqual({k: m.get(k) for k in expect}, expect)
        self.assertNotIn(('bcftools', 'refchanged', 'all', 'recall'), m)
        self.assertNotIn(('bcftools', 'refchanged', 'all', 'native_only_none'), m)


if __name__ == '__main__':
    unittest.main()
