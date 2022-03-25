import unittest
from pyHLAMSA import Genemsa
from tempfile import NamedTemporaryFile, mkstemp
from Bio import SeqIO, AlignIO


class TestMsaMainFunction(unittest.TestCase):
    """
    This test the core functionality
    """
    # TODO
    # format_alignment
    # format_alignment_diff
    # format_alignment_from_center
    # format_variantion_base
    # align

    def setUp(self):
        self.msa = Genemsa("yourname", seq_type="gen")
        alleles = {
            'a0': "CCATT-|GGT--GTCGGGT|TTC|C|AG",
            'a1': "CCACTG|GGT--ATCGGGT|TTC|C|AG",
            'c2': "CAATTG|GGT--GTCGGGT|---|A|AG",
        }
        self.msa.blocks = list(map(len, alleles['a0'].split("|")))
        self.msa.labels = [
            ["five_prime_UTR", "5UTR"],
            ["exon", "exon1"],
            ["intron", "intron1"],
            ["exon", "exon2"],
            ["three_prime_UTR", "3UTR"],
        ]
        self.msa.alleles = {k: v.replace("|", "") for k, v in alleles.items()}
        self.input_allele = alleles

    def test_shrink(self):
        # before shrink
        seq = self.input_allele['a1']
        self.assertEqual(self.msa.get_length(), len(seq.replace("|", "")))
        newmsa = self.msa.shrink()

        # after shrink
        # original msa not chanaged
        self.assertEqual(self.msa.blocks, list(map(len, seq.split("|"))))
        self.assertEqual(self.msa.get_length(), len(seq.replace("|", "")))
        # newmsa
        #    CCATT-|GGT--GTCGGGT|TTC|CAG
        # -> CCATT-|GGTGTCGGGT|TTC|CAG
        self.assertEqual(newmsa.get_length(), len(seq.replace("|", "").replace("-", "")))
        self.assertEqual(newmsa.blocks, list(map(len, seq.replace("-", "").split("|"))))

        # shrink c2 only: 0-length blocks because of gap
        c2_seq = self.input_allele['c2']
        c2 = self.msa.select_allele(["c2"]).shrink()
        self.assertEqual(c2.get_length(), len(c2_seq.replace("|", "").replace("-", "")))
        self.assertEqual(len(c2.blocks), c2_seq.count("|") + 1)
        self.assertEqual(c2.blocks[2], 0)

    def test_consensus(self):
        alleles = {k: v.replace("|", "") for k, v in self.input_allele.items()}

        # include_gap = True
        seq = self.msa.get_consensus(include_gap=True)
        self.assertTrue(len(seq), self.msa.get_length())
        for i, bases in enumerate(zip(*alleles.values())):
            self.assertTrue(list(bases).count(seq[i]) >= 0.5)

        # include_gap=False + gap exist
        seq = self.msa.get_consensus(include_gap=False)
        self.assertTrue(len(seq), self.msa.get_length())
        for i, bases in enumerate(zip(*alleles.values())):
            if not all(b == '-' for b in bases):
                self.assertTrue(bases.count(seq[i]) >= 0.5)
            else:
                self.assertTrue(seq[i] == "A")  # all gap => A, very special case
            self.assertTrue(seq[i] != "-")

        # include_gap=True + remove gap
        alleles = {k: v[:9] + v[11:] for k, v in alleles.items()}
        seq = self.msa.shrink().get_consensus(include_gap=False)
        for i, bases in enumerate(zip(*alleles.values())):
            self.assertTrue(bases.count(seq[i]) >= 0.5)
            self.assertTrue(seq[i] != "-")

    def test_save_fasta(self):
        fname = mkstemp()[1]

        # save fasta with gap
        n = 0
        self.msa.save_fasta(fname, gap=True)
        for seq in SeqIO.parse(fname, "fasta"):
            n += 1
            self.assertTrue(seq.id in self.input_allele)
            self.assertEqual(str(seq.seq), self.input_allele[seq.id].replace("|", ""))
        self.assertEqual(n, len(self.input_allele))

        # read gapped read by Biopython MSA
        msa = AlignIO.read(fname, "fasta")
        n = 0
        for seq in msa:
            n += 1
            self.assertTrue(seq.id in self.input_allele)
            self.assertEqual(str(seq.seq), self.input_allele[seq.id].replace("|", ""))
        self.assertEqual(n, len(self.input_allele))

        # no gap
        n = 0
        self.msa.save_fasta(fname, gap=False)
        for seq in SeqIO.parse(fname, "fasta"):
            n += 1
            self.assertTrue(seq.id in self.input_allele)
            self.assertEqual(str(seq.seq), self.input_allele[seq.id].replace("|", "").replace("-", ""))
        self.assertEqual(n, len(self.input_allele))

    def test_load_save_msa(self):
        fname1 = mkstemp()[1]
        fname2 = mkstemp()[1]
        self.msa.save_msa(fname1, fname2)
        newmsa = Genemsa.load_msa(fname1, fname2)

        # check for same msa
        self.assertEqual(len(newmsa), len(self.input_allele))
        self.assertEqual(newmsa.blocks, self.msa.blocks)
        self.assertEqual(newmsa.labels, self.msa.labels)
        for name in newmsa.get_sequence_names():
            self.assertEqual(newmsa.get(name), self.input_allele[name].replace("|", ""))

    def test_allele_selection(self):
        # select by list
        self.assertEqual(len(self.msa.select_allele(["a1", "a0"]).alleles), 2)
        with self.assertRaises(KeyError):
            self.msa.select_allele(["a1", "a2"])
        self.assertEqual(len(self.msa.select_allele([]).alleles), 0)

        # select by regex
        self.assertEqual(len(self.msa.select_allele(r"a.*").alleles), 2)
        self.assertEqual(len(self.msa.select_allele(r"b.*").alleles), 0)

        # select the sequence
        self.assertEqual(self.msa.get("a1"), self.input_allele["a1"].replace("|", ""))
        with self.assertRaises(KeyError):
            self.msa.get("a3")

        # change doesn't effect previous object
        self.msa.select_allele(r"a.*").alleles['a1'] = "123"
        self.assertNotEqual(self.msa.get('a1'), "123")
        self.msa.select_allele(r"a.*").blocks[1] = 123
        self.assertNotEqual(self.msa.blocks[1], 123)

    def test_add_seq(self):
        seq = self.input_allele['a1'].replace("|", "").replace("G", "A")
        n = len(self.input_allele)

        # add new sequence
        self.msa.add("test", seq)
        self.assertEqual(len(self.msa), n + 1)
        self.msa.add("test2", seq)
        self.assertEqual(len(self.msa), n + 2)

        # same name error
        with self.assertRaises(ValueError):
            self.msa.add("test", seq)
        self.assertEqual(len(self.msa), n + 2)

        #  different length
        with self.assertRaises(ValueError):
            self.msa.add("test1", seq + "C")
        self.assertEqual(len(self.msa), n + 2)

        # remove
        self.msa.remove(["test", "test2"])
        self.assertEqual(len(self.msa), n)

        # fail to merge msa (extend)
        newmsa = self.msa.copy()
        with self.assertRaises(ValueError):
            self.msa.extend(newmsa)
        with self.assertRaises(ValueError):
            newmsa.alleles["a1"] = "C"
            self.msa.extend(newmsa)

        # extend: success
        self.assertEqual(len(self.msa), n)
        newmsa.alleles = {"test3": self.input_allele["a1"].replace("|", "")}
        self.msa.extend(newmsa)
        self.assertEqual(len(self.msa), n + 1)

        # cleanup
        self.msa.remove("test3")
        self.assertEqual(len(self.msa), n)

    def test_reverse(self):
        newmsa = self.msa.reverse_complement()
        newmsa = newmsa.reverse_complement()
        for name in newmsa.get_sequence_names():
            self.assertEqual(newmsa.get(name), self.input_allele[name].replace("|", ""))
        self.assertEqual(self.msa.blocks, newmsa.blocks)

        # modify new msa wll nott change orignal msa
        newmsa.alleles['a1'] = "123"
        self.assertEqual(self.msa.get("a1"), self.input_allele["a1"].replace("|", ""))
        newmsa.blocks[1] = 123
        self.assertNotEqual(self.msa.blocks, newmsa)

    def test_basic_attr(self):
        n = len(self.input_allele)
        self.assertEqual(self.msa.get_sequence_names(), list(self.input_allele.keys()))
        # shape information checking
        self.assertEqual(self.msa.get_sequence_num(), n)
        self.assertEqual(len(self.msa), n)
        self.assertEqual(self.msa.get_length(), len(self.input_allele['a1'].replace("|", "")))
        self.assertEqual(self.msa.size(), (n, len(self.input_allele['a1'].replace("|", ""))))

        # test deepcopy
        newmsa = self.msa.copy()
        newmsa.alleles['a1'] = "123"
        self.assertEqual(self.msa.get("a1"), self.input_allele["a1"].replace("|", ""))
        newmsa.blocks[1] = 123
        self.assertNotEqual(self.msa.blocks, newmsa)

        # test Bio.MultipleSeqAlignment
        newmsa = Genemsa.from_MultipleSeqAlignment(self.msa.to_MultipleSeqAlignment())
        for name in newmsa.get_sequence_names():
            self.assertEqual(newmsa.get(name), self.input_allele[name].replace("|", ""))

    def test_select_chunk(self):
        chunk_num = self.input_allele['a1'].count("|") + 1
        self.assertEqual(self.msa.get_sequence_names(), list(self.input_allele.keys()))

        # check: out of range
        with self.assertRaises(IndexError):
            self.msa.select_exon([0])
        with self.assertRaises(IndexError):
            self.msa.select_exon([chunk_num // 2 + 1])
        with self.assertRaises(IndexError):
            self.msa.select_chunk([-2])
        with self.assertRaises(IndexError):
            self.msa.select_chunk([chunk_num])

        # chunk vs exon
        exon = self.msa.select_exon([1, 2])  # exon index start from 1
        exon_1 = self.msa.select_chunk([1, 3])  # chunk index start from 0
        self.assertEqual(len(exon), len(self.input_allele))
        self.assertEqual(exon.get_length(), exon_1.get_length())

        # check it's exon
        for label in exon.labels:
            self.assertEqual(label[0], "exon")
        for name in exon.get_sequence_names():
            self.assertEqual(exon.get(name), "".join(self.input_allele[name].split("|")[1::2]))

        # check select all exon
        exon_2 = self.msa.select_exon()
        self.assertEqual(exon.get_length(), exon_2.get_length())
        # check select last intron
        utr = self.msa.select_chunk([-1])
        for name in utr.get_sequence_names():
            self.assertEqual(utr.get(name), self.input_allele[name].split("|")[-1])

        # TODO: set type nuc or intron is not odd to test err

    def test_select_base(self):
        # use list to select base
        ind = [1, 2, 4]
        newmsa = self.msa[ind]
        for name in newmsa.get_sequence_names():
            for ind_i, i in enumerate(ind):
                self.assertEqual(newmsa.get(name)[ind_i], self.input_allele[name].replace("|", "")[i])

        # use slice to range of bases
        ind = slice(3, 10)
        newmsa = self.msa[ind]
        for name in newmsa.get_sequence_names():
            self.assertEqual(newmsa.get(name), self.input_allele[name].replace("|", "")[ind])

    def test_variantion(self):
        # select bases has variation
        bases = self.msa.get_variantion_base()
        newmsa = self.msa[bases]
        for i in zip(*newmsa.alleles.values()):
            self.assertTrue(len(set(i)) > 1)

        # bases with same sequence
        no_var_bases = sorted(set(range(self.msa.get_length())) - set(bases))
        newmsa = self.msa[no_var_bases]
        self.assertEqual(len(set(newmsa.alleles.values())), 1)

    def test_cigar(self):
        # I didn't not test save_bam
        # Match = 0
        # Insert = 1
        # Delete = 2
        # Mismatch = 8
        self.assertEqual(self.msa._calculate_cigar("AAA", "C--"), [[8, 1], [2, 2]])
        self.assertEqual(self.msa._calculate_cigar("AAA", "--C"), [[2, 2], [8, 1]])
        self.assertEqual(self.msa._calculate_cigar("AAA", "A-C"), [[0, 1], [2, 1], [8, 1]])
        self.assertEqual(self.msa._calculate_cigar("A-A", "ACC"), [[0, 1], [1, 1], [8, 1]])

    def test_gff(self):
        # because all sequences has same exon intron position
        # and already test it in load_msa and save_msa
        pass

    def test_split(self):
        msa_blocks = self.msa.split()
        self.assertEqual(self.input_allele['a1'].count("|") + 1, len(msa_blocks))
        self.assertEqual([i.get_length() for i in msa_blocks], list(map(len, self.input_allele['a1'].split("|"))))

        for newmsa in msa_blocks:
            self.assertEqual(len(newmsa.blocks), 1)
            self.assertEqual(len(newmsa.labels), 1)
            self.assertEqual(len(newmsa), len(self.input_allele))

    def test_concat(self):
        s = self.msa.split()
        newmsa = s[0]
        for i in s[1:]:
            newmsa += i
        print(newmsa)
        self.assertEqual(newmsa.blocks, self.msa.blocks)
        self.assertEqual(newmsa.labels, self.msa.labels)


class TestMsaExonOnly(unittest.TestCase):
    """
    This deal with the case where MSA contains exon-only sequences
    """
    def test_partical_seq(self):
        # setup
        msa = Genemsa("yourname", seq_type="gen")
        alleles = {
            'a0': "CCATT-|GGT--GTCGGGT|TTC|C|AG",
            'a1': "CCACTG|GGT--ATCGGGT|TTC|C|AG",
            'c2': "CAATTG|GGT--GTCGGGT|---|A|AG",
            'e0': "EEEEEE|GGT--ATCGGGT|EEE|C|AG",  # edit from a1
            'e1': "EEEEEE|GGT--ATCGGGT|EEE|E|AG",  # edit from a1
            'e2': "CCACTG|GGT--ATCGGGT|ETC|C|AG",  # edit from a1
        }
        msa.blocks = list(map(len, alleles['a0'].split("|")))
        msa.labels = [
            ["five_prime_UTR", "5UTR"],
            ["exon", "exon1"],
            ["intron", "intron1"],
            ["exon", "exon2"],
            ["three_prime_UTR", "3UTR"],
        ]
        msa.alleles = {k: v.replace("|", "") for k, v in alleles.items()}

        # Check E in sequence
        self.assertEqual(msa.select_complete().get_sequence_names(), ["a0", "a1", "c2"])
        self.assertEqual(msa.select_incomplete().get_sequence_names(), ["e0", "e1", "e2"])
        self.assertEqual(msa.select_exon().select_incomplete().get_sequence_names(), ["e1"])

        # check fill the E
        partical_seq_name = msa.select_incomplete().get_sequence_names()
        newmsa = msa.fill_incomplete('a1')
        for name in partical_seq_name:
            self.assertEqual(msa.get(name), msa.get('a1'))

    def test_merge_exon(self):
        # setup
        alleles = {
            'a0': "CCATT-|GGT--GTCGGGT|TTC|C|AG",
            'a1': "CCACTG|GGT--ATCGGGT|TTC|C|AG",
            'c2': "CAATTG|GGT--GTCGGGT|---|A|AG",
        }
        msa = Genemsa("yourname", seq_type="gen")
        msa.blocks = list(map(len, alleles['a0'].split("|")))
        msa.labels = [
            ["five_prime_UTR", "5UTR"],
            ["exon", "exon1"],
            ["intron", "intron1"],
            ["exon", "exon2"],
            ["three_prime_UTR", "3UTR"],
        ]
        msa.alleles = {k: v.replace("|", "") for k, v in alleles.items()}

        # selecct exon part and rename: add "e" before the name
        msa_nuc = msa.select_exon()
        msa_nuc.alleles.update({"e" + name: seq for name, seq in msa_nuc.alleles.items()})

        # main
        msa_merged = msa.merge_exon(msa_nuc)

        # test exon is same
        msa_merged_exon = msa_merged.select_incomplete().select_exon()
        self.assertEqual(len(msa_merged_exon), len(alleles))
        for name in msa_merged_exon.get_sequence_names():
            self.assertEqual(msa_merged_exon.get(name), msa_nuc.get(name))

        # test intron is E
        msa_merged_intron = msa_merged.select_incomplete().select_chunk(list(range(0, len(msa.blocks), 2)))
        for name in msa_merged_intron.get_sequence_names():
            self.assertEqual(set(msa_merged_intron.get(name)), set("E"))


if __name__ == '__main__':
    unittest.main()
