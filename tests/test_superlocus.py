import unittest
from collections import namedtuple
import fixtures
import vgraph.match as match
import vgraph.norm as norm

MockSample = namedtuple('MockSample', ['allele_indices', 'phased'])
MockRecord = namedtuple('MockRecord', ['contig', 'start', 'stop', 'ref', 'alleles', 'samples'])


class TestSuperlocusEquality(unittest.TestCase):

    def _make_loci(self, ref, vars):
        recs = [
            MockRecord(
                contig=contig,
                start=start,
                stop=start+len(alleles[0]),
                ref=ref,
                alleles=alleles,
                samples={'sample': MockSample(allele_indices=allele_indices, phased=False)})
            for contig, start, alleles, allele_indices in vars
            ]
        return [norm.NormalizedLocus(i, record, ref, name='sample') for i, record in enumerate(recs)]

    def test_snp(self):
        """Test the trivial matcher for minimal functionality"""
        ref = fixtures.FASTA_1

        sl1 = self._make_loci(ref, [('1', 2, ('T', 'A',), (0, 1))])
        sl2 = self._make_loci(ref, [('1', 2, ('T', 'A',), (0, 1))])

        matched, status = match.superlocus_equal(ref, None, None, sl1, sl2)
        self.assertTrue(matched)

    def test_complex1(self):
        """A fairly simple case that uses the haplotype matcher"""
        ref = fixtures.FASTA_1

        sl1 = self._make_loci(ref, [('1', 2, ('TGA', 'AGC',), (0, 1))])
        sl2 = self._make_loci(ref, [('1', 2, ('T', 'A',), (0, 1)),
                                    ('1', 4, ('A', 'C',), (0, 1))])
        matched, status = match.superlocus_equal(ref, None, None, sl1, sl2)
        self.assertTrue(matched)


    # def test_complex_dels(self):
    #     """A handful of deletions requiring normalization & complex matching"""
    #     ref = fixtures.FASTA_1
    #
    #     sl1 = self._make_loci(ref, [('1', 12, ('AGAG', '',), (0, 1))])
    #     sl2 = self._make_loci(ref, [('1', 9, ('A', '',), (0, 1)),
    #                                 ('1', 13, ('G', '',), (0, 1)),
    #                                 ('1', 16, ('AG', '',), (0, 1))])
    #     matched, status = match.superlocus_equal(ref, None, None, sl1, sl2)
    #     self.assertTrue(matched)

    def test_zygosity1(self):
        ref = fixtures.FASTA_1
        sl1 = self._make_loci(ref, [('1', 2, ('TGA', 'AGC',), (1, 1))])
        sl2 = self._make_loci(ref, [('1', 2, ('T', 'A',), (1, 1)),
                                    ('1', 4, ('A', 'C',), (0, 1))])
        matched, status = match.superlocus_equal(ref, None, None, sl1, sl2)

        self.assertFalse(matched)
        self.assertTrue(status == 'Z')

    def test_zygosity2(self):
        ref = fixtures.FASTA_1
        sl1 = self._make_loci(ref, [('1', 2, ('TGA', 'AGC',), (0, 1))])
        sl2 = self._make_loci(ref, [('1', 2, ('T', 'A',), (1, 1)),
                                    ('1', 4, ('A', 'C',), (1, 1))])
        matched, status = match.superlocus_equal(ref, None, None, sl1, sl2)

        self.assertFalse(matched)
        self.assertTrue(status == 'Z')

    def test_zygosity3(self):
        ref = fixtures.FASTA_1
        sl1 = self._make_loci(ref, [('1', 2, ('TGA', 'AGC',), (0, 1))])
        sl2 = self._make_loci(ref, [('1', 2, ('T', 'A',), (1, 1)),
                                    ('1', 4, ('A', 'C',), (0, 1))])
        matched, status = match.superlocus_equal(ref, None, None, sl1, sl2)

        self.assertFalse(matched)
        self.assertTrue(status == 'Z')


if __name__=="__main__":
    unittest.main()