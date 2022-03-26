import os
import unittest
from pyHLAMSA import Genemsa, Readmsa
from tempfile import NamedTemporaryFile, mkstemp
from Bio import SeqIO, AlignIO


class TestMsaReadFromDB(unittest.TestCase):
    """
    Test the correctness of reading data provided by Database
    """
    def setUp(self):
        # os.system("gunzip -k ./tests/hla_A01.dat.gz")
        # os.system("gunzip -k ./tests/A01_gen.txt.gz")
        # os.system("gunzip -k ./tests/A01_gen.msf.gz")
        pass

    def test_dat(self):
        data = Readmsa.read_dat_block("./tests/hla_A01.dat")
        self.assertEqual(data['A*01:01:01:01'], [
            ['UTR', 1, 300],
            ['exon1', 301, 373],
            ['intron1', 374, 503],
            ['exon2', 504, 773],
            ['intron2', 774, 1014],
            ['exon3', 1015, 1290],
            ['intron3', 1291, 1869],
            ['exon4', 1870, 2145],
            ['intron4', 2146, 2247],
            ['exon5', 2248, 2364],
            ['intron5', 2365, 2806],
            ['exon6', 2807, 2839],
            ['intron6', 2840, 2981],
            ['exon7', 2982, 3029],
            ['intron7', 3030, 3198],
            ['exon8', 3199, 3203],
            ['UTR', 3204, 3503]
        ])
        self.assertEqual(data['A*01:01:64'], [
            ['exon2', 1, 270],
            ['exon3', 271, 546]
        ])

    def test_msa_file(self):
        # this function will run parse_alignment also
        msa = Readmsa.from_alignment_file(f"tests/A01_gen.txt")
        msa.gene_name = "A"
        msa.seq_type = "gen"
        msa._assume_label()
        self.assertEqual(len(msa.blocks), 17)
        self.assertEqual(msa.get_length(), 3834)
        self.assertEqual(len(msa), 301)

    def test_msf_file(self):
        msa = Readmsa.from_MSF_file("tests/A01_gen.msf")
        self.assertEqual(len(msa.blocks), 1)
        self.assertEqual(msa.get_length(), 3890)
        self.assertEqual(len(msa), 301)
        # TODO: test merge_dat

    def test_same_in_two_method(self):
        # msf
        msa = Readmsa.from_MSF_file("tests/A01_gen.msf")
        # msa
        msa_ali = Readmsa.from_alignment_file(f"tests/A01_gen.txt")
        for name in msa.get_sequence_names():
            self.assertEqual(msa.get(name).replace('-', ''),
                             msa_ali.get(name).replace('-', ''))
        self.assertEqual(sorted(msa.get_sequence_names()), sorted(msa_ali.get_sequence_names))

    def test_msf_with_dat(self):
        # main test
        dat = Readmsa.read_dat_block("tests/hla_A01.dat")
        msa = Readmsa.from_MSF_file("tests/A01_gen.msf")
        msa1 = Readmsa.apply_dat_info_on_msa(msa, dat)
        self.assertEqual(len(msa1.blocks), 17)
        self.assertEqual(msa1.get_length(), 3834)
        self.assertEqual(len(msa1), 301)

        # compare result
        msa2 = Readmsa.from_alignment_file(f"tests/A01_gen.txt")
        self.assertEqual(sorted(msa1.get_sequence_names()), sorted(msa2.get_sequence_names))
        for m1, m2 in zip(msa1.split(), msa2.split()):
            for name in msa1.get_sequence_names():
                self.assertEqual(msa1.get(name).replace('-', ''),
                                 msa2.get(name).replace('-', ''))


if __name__ == '__main__':
    unittest.main()
