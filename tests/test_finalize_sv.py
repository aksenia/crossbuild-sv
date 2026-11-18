"""finalize_sv.py: END rules, chain QC reasons, REF check, end-to-end tool mode."""
import os
import subprocess
import sys
import tempfile
import unittest

import pysam

from liftover import finalize_sv as fs

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# chr1 source -> chr1 target, + strand, three blocks (0-based half-open):
#   src   0-300 -> tgt 100-400   then a 50 bp gap in the target only
#   src 300-500 -> tgt 450-650   then a 100 bp gap in the source only (500-600)
#   src 600-1000 -> tgt 650-1050
# chr2 source -> chr2 target on the - strand (inverted), one 500 bp block.
# chr3: a + chain with a 300 bp gap (300-600) holding a nested inverted chain,
# as UCSC chains represent an inversion.
CHAIN = """\
chain 1000 chr1 1000 + 0 1000 chr1 2000 + 100 1050 1
300 0 50
200 100 0
400

chain 500 chr2 1000 + 0 500 chr2 1000 - 0 500 2
500

chain 700 chr3 1000 + 0 1000 chr3 1000 + 0 1000 3
300 300 300
400

chain 300 chr3 1000 + 300 600 chr3 1000 - 400 700 4
300
"""
SEQ = 'ACGT' * 500   # hg38 chr1, 2000 bp; chr2 is its first 1000 bp


def body(path):
    with open(path) as fh:
        return [l.split('\t') for l in fh if not l.startswith('#')]


def parts(chrom, pos, alt, ref='A', vid='A'):
    return [chrom, str(pos), vid, ref, alt, '.', 'PASS', '.']


class FinalizeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        d = cls.tmp.name
        cls.chain_path = os.path.join(d, 'tiny.chain')
        with open(cls.chain_path, 'w') as fh:
            fh.write(CHAIN)
        cls.fasta_path = os.path.join(d, 'hg38.fa')
        with open(cls.fasta_path, 'w') as fh:
            fh.write('>chr1\n' + SEQ + '\n>chr2\n' + SEQ[:1000] + '\n')
        pysam.faidx(cls.fasta_path)
        cls.chains = fs.ChainIndex(cls.chain_path)

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def end(self, p, info, src=None):
        source = {'A': fs.SourceVariant(*src)} if src else {}
        return fs.target_end(p, info, source, self.chains)

    def reason(self, *args, **kw):
        with self.assertRaises(ValueError) as cm:
            self.end(*args, **kw)
        return str(cm.exception)

    def test_chain_mapping(self):
        m = self.chains.by_source['1'][0]
        self.assertEqual(m.map_point(101), 201)
        self.assertEqual(m.map_point(351), 501)    # second block, after target gap
        self.assertIsNone(m.map_point(550))        # source gap

    def test_sequence_end_from_ref(self):
        self.assertEqual(self.end(parts('chr1', 201, 'A', ref='ACGT'), {'SVTYPE': 'DEL'}), 204)

    def test_symbolic_ins_end_is_pos(self):
        self.assertEqual(self.end(parts('chr1', 201, '<INS>'), {'SVTYPE': 'INS'}), 201)

    def test_symbolic_interval_end_from_chain(self):
        self.assertEqual(self.end(parts('chr1', 201, '<DEL>'), {'SVTYPE': 'DEL'}, ('1', 101, 351)), 501)

    def test_end_not_in_chain(self):
        r = self.reason(parts('chr1', 201, '<DEL>'), {'SVTYPE': 'DEL'}, ('1', 101, 550))
        self.assertEqual(r, 'END_NOT_IN_CHAIN')

    def test_chain_does_not_match_tool_pos(self):
        r = self.reason(parts('chr1', 202, '<DEL>'), {'SVTYPE': 'DEL'}, ('1', 101, 351))
        self.assertEqual(r, 'CHAIN_DOES_NOT_MATCH_TOOL_POS')

    def test_inverted_chain(self):
        # src 101..201 maps to 900..800: POS must be the left end (800)
        self.assertEqual(self.end(parts('chr2', 800, '<DUP>'), {'SVTYPE': 'DUP'}, ('2', 101, 201)), 900)
        r = self.reason(parts('chr2', 900, '<DUP>'), {'SVTYPE': 'DUP'}, ('2', 101, 201))
        self.assertEqual(r, 'INVERTED_CHAIN')

    def test_non_base_alt(self):
        # transanno reverse-complements <DUP:TANDEM> on inverted chains
        r = self.reason(parts('chr2', 800, 'XXXXNTAXXXXX'), {'SVTYPE': 'DUP'}, ('2', 101, 201))
        self.assertEqual(r, 'NON_BASE_ALT')

    def test_native_mode_keeps_end(self):
        p = parts('chr1', 201, '<DEL>')
        self.assertEqual(fs.target_end(p, {'SVTYPE': 'DEL', 'END': '300'}, {}, None), 300)
        with self.assertRaises(ValueError):
            fs.target_end(p, {'SVTYPE': 'DEL'}, {}, None)

    def test_check_ref(self):
        fasta = pysam.FastaFile(self.fasta_path)
        fs.check_ref(fasta, 'chr1', 201, 'ACGT')
        for chrom, ref, reason in (('chr1', 'C', 'REF_MISMATCH'), ('chrZ', 'A', 'CONTIG_NOT_IN_HG38')):
            with self.assertRaises(ValueError) as cm:
                fs.check_ref(fasta, chrom, 201, ref)
            self.assertEqual(str(cm.exception), reason)

    def test_lift_strand(self):
        src = fs.SourceVariant('1', 101, 101)
        self.assertEqual(self.chains.strand(src, 'chr1', 201), '+')
        self.assertEqual(self.chains.strand(src, 'chr2', 201), '.')     # no chain to chr2
        # anchor base in the source gap (500-600): no block, but inside the chain span
        self.assertIsNone(self.chains.by_source['1'][0].map_point(550))
        self.assertEqual(self.chains.strand(fs.SourceVariant('1', 550, 550), 'chr1', 651), '+')
        self.assertEqual(self.chains.strand(fs.SourceVariant('2', 101, 201), 'chr2', 800), '-')
        # nested inversion: the + chain spans it but has no block there
        self.assertEqual(self.chains.strand(fs.SourceVariant('3', 401, 401), 'chr3', 500), '-')

    def test_oriented_ci(self):
        ci = {'CIPOS': '-10,2', 'CIEND': '-3,7'}
        self.assertEqual(fs.oriented_ci(ci, '+', 'DEL'), ci)
        self.assertEqual(fs.oriented_ci(ci, '-', 'DEL'), {'CIPOS': '-7,3', 'CIEND': '-2,10'})
        self.assertEqual(fs.oriented_ci({'CIPOS': '-10,2'}, '-', 'DEL'), {'CIPOS': '0,0', 'CIEND': '-2,10'})
        self.assertEqual(fs.oriented_ci({'CIPOS': '-10,2'}, '-', 'INS'), {'CIPOS': '-2,10'})

    def test_main_tool_mode(self):
        d = self.tmp.name
        header = '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
        source = os.path.join(d, 'source.vcf')
        with open(source, 'w') as fh:
            fh.write(header + '1\t101\tA\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=351\n')
        vcf = header + ('chr1\t201\tA\tA\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=999\n'
                        'chr1\t202\tB\tA\tAT\t.\tPASS\tSVTYPE=INS\n')    # hg38 base is C
        out, rej = os.path.join(d, 'out.vcf'), os.path.join(d, 'rej.vcf')
        subprocess.run([sys.executable, os.path.join(ROOT, 'liftover', 'finalize_sv.py'),
                        '--tool', 'crossmap', '--source', source, '--chain', self.chain_path,
                        '--fasta', self.fasta_path, '--output', out, '--reject', rej],
                       input=vcf, text=True, check=True, capture_output=True)
        kept, rejected = body(out), body(rej)
        self.assertEqual([p[2] for p in kept], ['A'])
        self.assertEqual(kept[0][7].strip(), 'SVTYPE=DEL;END=501;LIFT_STRAND=+;LIFTOVER_TOOL=crossmap')
        self.assertEqual(rejected[0][6], 'QC_FAILED')
        self.assertIn('QC_REJECT_REASON=REF_MISMATCH', rejected[0][7])

    def test_ci_from_source(self):
        # CrossMap lifts CIEND like a position; the source offsets are restored,
        # and on the inverted chain (chr2) negated and swapped
        d = self.tmp.name
        header = '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
        source = os.path.join(d, 'source_ci.vcf')
        with open(source, 'w') as fh:
            fh.write(header + '1\t101\tA\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=351;CIPOS=-5,5;CIEND=0,1\n'
                     '2\t101\tB\tN\t<DUP>\t.\tPASS\tSVTYPE=DUP;END=201;CIPOS=-10,2;CIEND=-3,7\n')
        vcf = header + ('chr1\t201\tA\tA\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=501;CIPOS=-5,5;CIEND=501,1\n'
                        'chr2\t800\tB\tT\t<DUP>\t.\tPASS\tSVTYPE=DUP;END=900;CIPOS=-10,2;CIEND=-3,7\n')
        out, rej = os.path.join(d, 'out_ci.vcf'), os.path.join(d, 'rej_ci.vcf')
        subprocess.run([sys.executable, os.path.join(ROOT, 'liftover', 'finalize_sv.py'),
                        '--tool', 'crossmap', '--source', source, '--chain', self.chain_path,
                        '--fasta', self.fasta_path, '--output', out, '--reject', rej],
                       input=vcf, text=True, check=True, capture_output=True)
        self.assertEqual([p[7].strip() for p in body(out)], [
            'SVTYPE=DEL;END=501;CIPOS=-5,5;CIEND=0,1;LIFT_STRAND=+;LIFTOVER_TOOL=crossmap',
            'SVTYPE=DUP;END=900;CIPOS=-7,3;CIEND=-2,10;LIFT_STRAND=-;LIFTOVER_TOOL=crossmap'])


if __name__ == '__main__':
    unittest.main()
