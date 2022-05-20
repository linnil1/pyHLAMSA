import os
import unittest
from pyhlamsa import Genemsa, msaio
from pyhlamsa.utils import dat as datop

from tempfile import NamedTemporaryFile, mkstemp
from Bio import SeqIO, AlignIO


class TestMsaReadFromDB(unittest.TestCase):
    """
    Test the correctness of reading data provided by Database
    """
    def setUp(self):
        os.system("gunzip -k ./tests/hla_A01.dat.gz")
        os.system("gunzip -k ./tests/A01_gen.txt.gz")
        os.system("gunzip -k ./tests/A01_gen.msf.gz")
        pass

    def test_dat(self):
        data = datop.read_dat_block("./tests/hla_A01.dat")
        self.assertEqual(data['A*01:01:01:01'], [
            {'name': 'UTR', 'start': 1, 'end': 300},
            {'name': 'exon1', 'start': 301, 'end': 373},
            {'name': 'intron1', 'start': 374, 'end': 503},
            {'name': 'exon2', 'start': 504, 'end': 773},
            {'name': 'intron2', 'start': 774, 'end': 1014},
            {'name': 'exon3', 'start': 1015, 'end': 1290},
            {'name': 'intron3', 'start': 1291, 'end': 1869},
            {'name': 'exon4', 'start': 1870, 'end': 2145},
            {'name': 'intron4', 'start': 2146, 'end': 2247},
            {'name': 'exon5', 'start': 2248, 'end': 2364},
            {'name': 'intron5', 'start': 2365, 'end': 2806},
            {'name': 'exon6', 'start': 2807, 'end': 2839},
            {'name': 'intron6', 'start': 2840, 'end': 2981},
            {'name': 'exon7', 'start': 2982, 'end': 3029},
            {'name': 'intron7', 'start': 3030, 'end': 3198},
            {'name': 'exon8', 'start': 3199, 'end': 3203},
            {'name': 'UTR', 'start': 3204, 'end': 3503}
        ])

        self.assertEqual(data['A*01:01:64'], [
            {'name': 'exon2', 'start': 1, 'end': 270},
            {'name': 'exon3', 'start': 271, 'end': 546},
        ])

    def test_msa_file(self):
        # this function will run parse_alignment also
        msa = msaio.read_alignment_txt(f"tests/A01_gen.txt")
        msa.gene_name = "A"
        msa.assume_label("gen")
        self.assertEqual(len(msa.blocks), 17)
        self.assertEqual(msa.get_length(), 3834)
        self.assertEqual(len(msa), 301)

    def test_msf_file(self):
        msa = msaio.read_msf_file("tests/A01_gen.msf")
        self.assertEqual(len(msa.blocks), 1)
        self.assertEqual(msa.get_length(), 3890)
        self.assertEqual(len(msa), 301)
        # TODO: test merge_dat

    def test_same_in_two_method(self):
        # msf
        msa = msaio.read_msf_file("tests/A01_gen.msf")
        # msa
        msa_ali = msaio.read_alignment_txt(f"tests/A01_gen.txt")
        for name in msa.get_sequence_names():
            self.assertEqual(msa.get(name).replace('-', ''),
                             msa_ali.get(name).replace('-', ''))
        self.assertEqual(sorted(msa.get_sequence_names()), sorted(msa_ali.get_sequence_names()))

    def test_msf_with_dat(self):
        # main test
        dat = datop.read_dat_block("tests/hla_A01.dat")
        msa = msaio.read_msf_file("tests/A01_gen.msf")
        msa.gene_name = "A"

        msa1 = datop.apply_dat_info_on_msa(msa, dat, seq_type="gen")
        self.assertEqual(len(msa1.blocks), 17)
        self.assertEqual(msa1.get_length(), 3890)

        # why A*01:11N be ignored
        # exon3:
        # in dat: 252 bases
        # exon            1015..1266
        # in msf: 276 bases
        # The sequence length in exon3 is not match
        self.assertEqual(len(msa1), 301 - 1)
        self.assertTrue("A*01:11N" not in msa1.get_sequence_names())

        # compare result
        msa2 = msaio.read_alignment_txt(f"tests/A01_gen.txt")
        self.assertEqual(sorted(set(msa1.get_sequence_names()) | set(["A*01:11N"])), sorted(msa2.get_sequence_names()))
        for m1, m2 in zip(msa1.split(), msa2.split()):
            for name in msa1.get_sequence_names():
                self.assertEqual(msa1.get(name).replace('-', ''),
                                 msa2.get(name).replace('-', ''))


class TestMsaHLA(unittest.TestCase):
    @unittest.skipIf(not os.path.exists("alignment_v3470"), "Tested on local")
    def test_hla_align(self):
        from pyhlamsa import HLAmsa, HLAmsaEX
        genes = ["A", "B", "C"]
        hla0 = HLAmsa(genes, filetype="gen", imgt_alignment_folder="alignment_v3470", version="3470")
        hla1 = HLAmsa(genes, imgt_alignment_folder="alignment_v3470")
        hla2 = HLAmsaEX(genes, imgt_folder="IMGT_v3470", version="3470")

        # test squence is correct
        for gene in genes:
            for name in set(hla0[gene].get_sequence_names()) & set(hla1[gene].get_sequence_names()):
                self.assertEqual(hla0[gene].get(name).replace("-", "").replace("E", ""),
                                 hla1[gene].get(name).replace("-", "").replace("E", ""))
            for name in set(hla1[gene].get_sequence_names()) & set(hla2[gene].get_sequence_names()):
                self.assertEqual(hla1[gene].get(name).replace("-", "").replace("E", ""),
                                 hla2[gene].get(name).replace("-", "").replace("E", ""))

        # test block squence is correct
        for gene in genes:
            self.assertEqual(len(hla0[gene].blocks), len(hla1[gene].blocks))
            self.assertEqual(len(hla0[gene].blocks), len(hla2[gene].blocks))
            for i in range(len(hla0[gene].blocks)):
                sub_hla0 = hla0[gene].select_block([i])
                sub_hla1 = hla1[gene].select_block([i])
                sub_hla2 = hla2[gene].select_block([i])
                for name in set(hla0[gene].get_sequence_names()) & set(hla1[gene].get_sequence_names()):
                    self.assertEqual(sub_hla0.get(name).replace("-", "").replace("E", ""),
                                     sub_hla1.get(name).replace("-", "").replace("E", ""))
                for name in set(hla1[gene].get_sequence_names()) & set(hla2[gene].get_sequence_names()):
                    self.assertEqual(sub_hla1.get(name).replace("-", "").replace("E", ""),
                                     sub_hla2.get(name).replace("-", "").replace("E", ""))

    @unittest.skipIf(not os.path.exists("KIR_v2100"), "Tested on local")
    def test_kir(self):
        from pyhlamsa import KIRmsa
        kir0 = KIRmsa(filetype=["gen"], ipd_folder="KIR_v2100", version="2100")
        kir1 = KIRmsa(ipd_folder="KIR_v2100")
        self.assertEqual(set(kir0.list_genes()), set(kir1.list_genes()))
        genes = set(kir0.list_genes())
        for gene in genes:
            for name in set(kir0[gene].get_sequence_names()) & set(kir1[gene].get_sequence_names()):
                self.assertEqual(kir0[gene].get(name).replace("-", "").replace("E", ""), kir1[gene].get(name).replace("-", "").replace("E", ""))

    @unittest.skipIf(not os.path.exists("pharmvar-5.1.10"), "Tested on local")
    def test_cyp(self):
        from pyhlamsa import CYPmsa
        cpy = CYPmsa(pharmvar_folder="pharmvar-5.1.10")
        # already lots of assert in CYPmsa._read_db_gene


if __name__ == '__main__':
    unittest.main()
