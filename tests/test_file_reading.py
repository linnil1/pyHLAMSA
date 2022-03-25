import unittest
from pyHLAMSA import Genemsa
from tempfile import NamedTemporaryFile, mkstemp
from Bio import SeqIO, AlignIO


class TestMsaReadFromDB(unittest.TestCase):
    """
    Test the correctness of reading data provided by Database
    """
    def test_dat(self):
        data = Genemsa.read_dat("./tests/hla.dat")
        self.assertEqual(data, {
            'A*01:01:01:01': [
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
            ],
            'A*01:01:64': [
                ['exon2', 1, 270],
                ['exon3', 271, 546]
            ]
        })

    def test_msa_file(self):
        # this function will run parse_alignment also
        msa = Genemsa("A", seq_type="gen")
        msa.read_alignment_file(f"tests/A_gen.txt")
        self.assertEqual(len(msa.blocks), 3)
        self.assertEqual(msa.get_length(), 400)
        self.assertEqual(len(msa), 301)

    def test_msf_file(self):
        msa = Genemsa("A", "gen")
        msa = msa.read_MSF_file("tests/A_gen.msf")
        self.assertEqual(len(msa.blocks), 1)
        self.assertEqual(msa.get_length(), 300)
        self.assertEqual(len(msa), 301)
        # TODO: test merge_dat

    def test_same_in_two_method(self):
        # msf
        msa = Genemsa("A", "gen")
        msa = msa.read_MSF_file("tests/A_gen.msf")

        # msa
        msa_ali = Genemsa("A", seq_type="gen")
        msa_ali.read_alignment_file(f"tests/A_gen.txt")
        for name in msa.get_sequence_names():
            seq1 = msa.get(name).replace('-', '')
            seq2 = msa_ali.get(name).replace('-', '')
            leng = min(len(seq1), len(seq2))
            self.assertEqual(seq1[:leng], seq2[:leng])


if __name__ == '__main__':
    unittest.main()
