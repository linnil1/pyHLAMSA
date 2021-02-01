import logging
import os
from glob import glob
import re
from pprint import pprint

from Bio.Align import MultipleSeqAlignment, PairwiseAligner
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO
import pysam

# setup logging
logger = logging.getLogger("pyIMGTHLA")
logger.setLevel(logging.DEBUG)
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)
logger.addHandler(ch)


class HLAmsa:
    def __init__(self, genes=[], filetype=["gen", "nuc"],
                 imgt_folder="alignments", version="3430"):
        self.imgt_folder = imgt_folder
        if not os.path.exists(self.imgt_folder):
            self.download(True, version=version)
        else:
            logger.info(f"IMGT ver={version} exists")
        assert os.path.exists(self.imgt_folder)

        if not genes:
            fs = glob(f"{self.imgt_folder}/*.txt")
            genes = set([f.split("/")[-1].split("_")[0] for f in fs])
            # TODO:
            # * E last exon not shown in nuc
            # * P has gen but no nuc exist?
            # * N the nuc has one more bp than in gen
            genes = genes - set(["ClassI", "DRB", "E", "P", "N"])
            genes = sorted(list(genes))

        logger.info(f"Read Gene {genes}")
        self.genes = {}
        for gene in genes:
            logger.info(f"Reading {gene}")
            self.genes[gene] = self.read_alignments(gene, filetype)
            logger.debug(f"Merged {self.genes[gene]}")

    def read_alignments(self, gene, filetype):
        """ read gene_filetype """
        if "gen" in filetype:
            msa_gen = Genemsa(gene, seq_type="gen")
            msa_gen.read_file(f"{self.imgt_folder}/{gene}_gen.txt")
            logger.debug(f"{msa_gen}")
        if "nuc" in filetype:
            msa_nuc = Genemsa(gene, seq_type="nuc")
            # Special Case: DRB* nuc are in DRB_nuc.txt
            if gene.startswith("DRB"):
                msa_nuc.read_file(f"{self.imgt_folder}/DRB_nuc.txt")
                logger.debug(f"DRB: {msa_nuc}")
                msa_nuc = msa_nuc.select_allele(gene + ".*")
                logger.debug(f"{msa_nuc}")
            else:
                msa_nuc.read_file(f"{self.imgt_folder}/{gene}_nuc.txt")
                logger.debug(f"{msa_nuc}")

        if "gen" in filetype and "nuc" in filetype:
            return msa_gen.merge_exon(msa_nuc)
        elif "gen" in filetype:
            return msa_gen
        elif "nuc" in filetype:
            return msa_nuc
        else:
            return None

    def download(self, download=True, version="3430"):
        """
        Download the IMGTHLA alignments folder to `imgt_folder`
        """
        if os.path.exists(self.imgt_folder):
            return
        # TODO: Auto find the latest version
        fname = f"Alignments_Rel_{version}"
        url = "ftp://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/" + fname + ".zip"
        logger.info(f"Download IMGT data from {url}")
        os.system(f"wget {url}")
        os.system(f"unzip {fname}.zip")
        os.system(f"mv alignments {self.imgt_folder}")


class Genemsa:
    def __init__(self, gene_name, seq_type=""):
        self.gene_name = gene_name
        self.alleles = {}
        self.blocks = []  # intron exon length
        self.seq_type = seq_type

    def __str__(self):
        return f"<{self.gene_name} {self.seq_type} "\
               f"alleles={len(self.alleles)} length={self.blocks}>"

    # some helper functions
    def calculate_position(self):
        pos = [0]
        for length in self.blocks:
            pos.append(pos[-1] + length)
        return pos

    def get_length(self):
        return sum(self.blocks)

    def get_first(self):
        return next(iter(self.alleles.items()))

    # reading functions
    def read_file(self, fname):
        """ read MSA format in IMGT alignments folder """
        alleles = self.parse_alignment(fname)
        for allele, seq in alleles.items():
            self.alleles[allele] = seq.replace("|", "")
            # calculate length of exons and introns
            if not self.blocks:
                self.blocks = [len(seq) for seq in seq.split("|")]

    def parse_alignment(self, fname):
        """ Read alignments file with IMGT-alignment format """
        # parse aligments
        alleles = {}
        ref_allele = ""
        for line in open(fname):
            line = line.strip()
            if not re.match(r"^\w+\*", line):
                continue

            match = re.findall(r"^(.*?) +(.*)", line)
            # assert match
            # emtpy seq (B_nuc.txt)
            if not match:
                continue
            allele, seq = match[0]

            if allele not in alleles:
                if not ref_allele:
                    ref_allele = allele
                alleles[allele] = ""

            alleles[allele] += seq

        # check sequences and replace
        rm_allele = []
        for allele in alleles:
            alleles[allele] = alleles[allele].replace(" ", "").replace("*", ".")
            if allele == ref_allele:
                alleles[allele] = alleles[allele].replace(".", "-")
                continue
            ref_seq = alleles[ref_allele]
            if len(alleles[allele]) != len(ref_seq):
                rm_allele.append(allele)
                continue
            # assert len(alleles[allele]) == len(ref_seq)

            seq = list(alleles[allele])
            for i in range(len(seq)):
                if seq[i] == "-":
                    seq[i] = ref_seq[i]

                if seq[i] == "|":
                    assert seq[i] == ref_seq[i]
                else:
                    assert seq[i] in "ATCG."
            alleles[allele] = ''.join(seq).replace(".", "-")

        for allele in rm_allele:
            logger.warning(f"Remove {allele} due to length in {fname}")
            del alleles[allele]
        return alleles

    # Position selection functions
    def select_exon(self, exon_index=[]):
        """
        Extract exon from gene

        exon_index: list
            The index of exon you want to extract
            Index start from 1
            e.g.
              1 for exon1
              2 for exon2
            Leave empty if you want all the exons
        """
        if self.seq_type != "gen":
            raise TypeError("Check this object is gen or not")
        if len(self.blocks) % 2 != 1:
            raise ValueError("introns + exon should be odd")

        # If not specific the index, extract all exons
        if not exon_index:
            exon_index = range(1, len(self.blocks), 2)
        else:
            if not (max(exon_index) <= len(self.blocks) // 2 and min(exon_index) > 0):
                raise ValueError("Check chunk index is correct")
            exon_index = [i * 2 - 1 for i in exon_index]

        new_msa = self.select_chunk(exon_index)
        new_msa.seq_type = "nuc"
        for allele, seq in new_msa.alleles.items():
            assert "E" not in seq
        return new_msa

    def select_chunk(self, index=[]):
        """
        Extract intron and exon from gene

        index: list
            The index of exon you want to extract
            Index start from 1
            Leave empty if you want all the exons
            e.g.
              0 for 5-UTR
              1 for exon1
              2 for intron1
              3 for exon2
              4 for 3-UTR(for two exons gene)
              -1 for 3-UTR(for all case)
        """
        if not index:
            return self

        if not (max(index) <= len(self.blocks) and min(index) >= -1):
            raise ValueError("Check chunk index is correct")

        # replace -1 (3-UTR)
        for i in range(len(index)):
            if index[i] == -1:
                index[i] = len(self.blocks) - 1

        # new a msa object
        new_msa = Genemsa(self.gene_name)
        for i in index:
            new_msa.blocks.append(self.blocks[i])

        # extract
        gen_pos = self.calculate_position()
        for allele, gen_seq in self.alleles.items():
            new_seq = "".join([gen_seq[gen_pos[i]:gen_pos[i + 1]] for i in index])
            new_msa.alleles[allele] = new_seq

        return new_msa

    def __getitem__(self, index=[]):
        """ select region in the sequences """
        if not index:
            return self

        # Extract specific region in alignment
        if isinstance(index, slice) or isinstance(index, int):
            new_msa = Genemsa(self.gene_name)
            new_msa.alleles = {allele: seq[index]
                               for allele, seq in self.alleles.items()}
            new_msa.blocks = [len(new_msa.get_first()[1])]
            return new_msa

        elif isinstance(index, tuple) or isinstance(index, list):
            new_msa = Genemsa(self.gene_name)
            new_msa.blocks = [len(index)]
            new_msa.alleles = {allele: ''.join([seq[i] for i in index])
                               for allele, seq in self.alleles.items()}
            return new_msa
        # Fail
        else:
            raise TypeError("Bad usage")

    # Allele selection functions
    def get(self, ref_allele):
        """ Extract specific allele in alignment """
        if ref_allele not in self.alleles:
            raise ValueError(f"{ref_allele} not found")
        return self.alleles[ref_allele]

    def select_allele(self, regex):
        """ Select allele by regex """
        new_msa = Genemsa(self.gene_name, self.seq_type)
        new_msa.blocks = self.blocks
        new_msa.alleles = {allele: seq for allele, seq in self.alleles.items()
                           if re.match(regex, allele)}
        return new_msa

    def select_complete(self):
        """ Select non gDNA sequences """
        new_msa = Genemsa(self.gene_name, self.seq_type)
        new_msa.blocks = self.blocks
        new_msa.alleles = {allele: seq for allele, seq in self.alleles.items()
                           if "E" not in seq}
        return new_msa

    def select_imcomplete(self):
        """ Select cDNA-only sequences """
        new_msa = Genemsa(self.gene_name, self.seq_type)
        new_msa.blocks = self.blocks
        new_msa.alleles = {allele: seq for allele, seq in self.alleles.items()
                           if "E" in seq}
        return new_msa

    # Sequence adding functions
    def add(self, name, seq):
        """ Add sequence into msa """
        if len(seq) != self.get_length():
            raise ValueError("Length not match to alignments")
        if not len(self.blocks):
            raise ValueError("MSA is empty")
        if name in self.alleles:
            raise ValueError(f"{name} already exist")

        self.alleles[name] = seq
        return self

    def extend(self, msa):
        """ Add msa into this msa """
        if self.blocks != msa.blocks:
            raise ValueError("Length is different")
        self.alleles.update(msa.alleles)
        return self

    # Functions deal with exon-only alleles
    def fill_imcomplete(self, ref_allele):
        """ Fill cDNA-only sequences with ref_allele sequences"""
        if ref_allele not in self.alleles:
            raise ValueError(f"{ref_allele} not found")

        ref_seq = self.alleles[ref_allele]
        for allele, seq in self.alleles.items():
            if "E" in seq:
                self.alleles[allele] = ''.join([seq[i] if seq[i] != "E" else ref_seq[i] for i in range(len(seq))])
        return self

    def merge_exon(self, msa_nuc):
        """ Merge nuc into gen """
        # merge when allele in both nuc and gen
        if not (self.seq_type == "gen" and msa_nuc.seq_type == "nuc"):
            raise TypeError("Should merge nuc into gen")
        if (self.gene_name != msa_nuc.gene_name or
                len(self.blocks) != len(msa_nuc.blocks) * 2 + 1):
            raise ValueError("Check object's name and chunks are correct")

        nuc_pos = msa_nuc.calculate_position()
        gen_pos = self.calculate_position()
        new_msa = Genemsa(self.gene_name, self.seq_type)
        ref_blocks = []

        # for each allele
        for allele, nuc_seq in msa_nuc.alleles.items():
            # check and merge into exons regions
            if allele in self.alleles:
                gen_seq = self.alleles[allele]
                new_seq = ""
                blocks = []
                for i in range(len(self.blocks)):
                    gen_seq_part = gen_seq[gen_pos[i]:gen_pos[i + 1]]
                    # intron and UTR
                    if i % 2 == 0:
                        new_seq += gen_seq_part
                        blocks.append(len(gen_seq_part))
                        continue

                    # exon
                    nuc_seq_part = nuc_seq[nuc_pos[i // 2]:nuc_pos[i // 2 + 1]]
                    if nuc_seq_part.replace("-", "") != gen_seq_part.replace("-", ""):
                        # Special Case: ignore exon5 in DQB1
                        if not (i == 9 and allele.startswith("DQB1")):
                            # print(nuc_seq_part)
                            # print(gen_seq_part)
                            logger.warning(f"Remove {allele} due to gen is different from other gen after insert nuc")
                            break
                    new_seq += nuc_seq_part
                    blocks.append(len(nuc_seq_part))
                else:  # no break triggered
                    # check and add to new msa
                    if not ref_blocks:
                        ref_blocks = blocks
                    else:
                        assert ref_blocks == blocks
                    new_msa.alleles[allele] = new_seq

            # Add exon only allele
            else:
                new_seq = ""
                for i in range(len(self.blocks)):
                    if i % 2 == 0:
                        new_seq += 'E' * self.blocks[i]
                    else:
                        new_seq += nuc_seq[nuc_pos[i // 2]:nuc_pos[i // 2 + 1]]
                new_msa.alleles[allele] = new_seq

        # Check
        non_insert_allele = set(self.alleles.keys()) - set(msa_nuc.alleles.keys())
        if len(non_insert_allele):
            if ref_blocks == self.blocks:
                for allele in non_insert_allele:
                    new_msa.alleles[allele] = self.alleles[allele]
            else:
                logger.warning(f"Remove {non_insert_allele} due to gen is different from other gen after insert nuc")

        new_msa.blocks = ref_blocks
        return new_msa

    # Format function
    def format_alignment_diff(self, ref_allele=""):
        """ Print all alleles diff from ref_allele """
        if not len(self.alleles):
            raise ValueError("MSA is empty")
        if not ref_allele:
            ref_allele = self.get_first()[0]
        if ref_allele not in self.alleles:
            raise ValueError(f"{ref_allele} not found")
        ref_seq = self.alleles[ref_allele]

        # use new msa object to save sequences
        new_msa = Genemsa(self.gene_name)
        new_msa.blocks = self.blocks
        new_msa.alleles = {ref_allele: ref_seq}
        for allele, seq in self.alleles.items():
            if allele == ref_allele:
                continue
            new_seq = ""
            for i in range(len(seq)):
                if seq[i] == ref_seq[i]:
                    if ref_seq[i] == '-':
                        new_seq += "."
                    else:
                        new_seq += "-"
                elif seq[i] == "-":
                    new_seq += "*"
                else:
                    new_seq += seq[i]
            new_msa.alleles[allele] = new_seq
        new_msa.alleles[ref_allele] = new_msa.alleles[ref_allele].replace("-", ".")
        return new_msa.format_alignment()

    def format_alignment(self, wrap=100):
        """ Format to xx_gen.txt """
        if not self.blocks or not self.alleles:
            raise ValueError("MSA is empty")
        output_str = ""
        bid = 0
        block_pos = self.calculate_position()
        block_pos, seq_length = block_pos[1:-1], block_pos[-1]

        for pos in range(0, seq_length, wrap):
            # Find the insertion point
            seq_part_len = min(seq_length, pos + wrap) - pos
            split_data = [(seq_part_len, 3, "\n")]  # pos, priority, char
            split_data += [(i, 2, " ") for i in range(10, seq_part_len, 10)]
            while bid < len(block_pos) and pos < block_pos[bid] <= pos + seq_part_len:
                split_data.append((block_pos[bid] - pos, 1, "|"))
                bid += 1
            split_data = sorted(split_data)

            # Header per wrapping
            output_str += f" {'gDNA':<18} {pos}\n"
            output_str += " " * 20 + "|\n"

            # alleles
            for allele, seq in self.alleles.items():
                output_str += f" {allele:18} "
                cur_pos = pos
                for j in split_data:
                    if j[0] > len(seq):
                        break
                    output_str += seq[cur_pos:pos + j[0]]
                    output_str += j[2]
                    cur_pos = pos + j[0]
            output_str += "\n\n"
        return output_str

    @classmethod
    def from_MultipleSeqAlignment(cls, bio_msa):
        """ Transfer MultipleSeqAlignment in biopython to this object """
        new_msa = Genemsa("", "")
        new_msa.blocks = [bio_msa.get_alignment_length()]
        for seq in bio_msa:
            new_msa.alleles[seq.id] = str(seq.seq)
            if len(seq.seq) != new_msa.blocks[0]:
                raise ValueError("Length is different inside MultipleSeqAlignment")
        return new_msa

    def to_MultipleSeqAlignment(self):
        """ Transfer this object to MultipleSeqAlignment in biopython """
        return MultipleSeqAlignment(self.to_fasta(gap=True))

    def to_fasta(self, gap=True):
        """ Transfer to list of SeqRecord """
        if gap:
            return [SeqRecord(Seq(seq.replace("E", "-")), id=allele, description="")
                    for allele, seq in self.alleles.items()]
        else:
            return [SeqRecord(Seq(seq.replace("E", "").replace("-", "")), id=allele, description="")
                    for allele, seq in self.alleles.items()]

    # Consensus function
    def calculate_frequency(self):
        """ Calculate ATCG and gap frequency in msa """
        freqs = []
        for i in zip(*self.alleles.values()):
            freqs.append([
                i.count("A"),
                i.count("T"),
                i.count("C"),
                i.count("G"),
                i.count("-")])
        return freqs

    def get_consensus(self, include_gap=False):
        """ Use max frequency one as consensus """
        freqs = self.calculate_frequency()
        if not include_gap:
            max_ind = [max(range(4), key=lambda i: f[i]) for f in freqs]
        else:
            max_ind = [max(range(5), key=lambda i: f[i]) for f in freqs]
        seq = ["ATCG-"[i] for i in max_ind]
        return ''.join(seq)

    def shrink(self):
        """ Remove empty base-pair across all allele """
        # index to delete
        freqs = self.calculate_frequency()
        masks = [f[4] != sum(f) for f in freqs]

        # recalcuate blocks
        new_msa = Genemsa(self.gene_name, self.seq_type)
        gen_pos = self.calculate_position()
        for i in range(len(self.blocks)):
            new_msa.blocks.append(sum(masks[gen_pos[i]:gen_pos[i+1]]))

        # remove bp in allele
        for allele, seq in self.alleles.items():
            new_msa.alleles[allele] = "".join([seq[i] for i in range(len(seq)) if masks[i]])
            assert len(new_msa.alleles[allele]) == new_msa.get_length()

        return new_msa

    # Align sequences on allele
    def align(self, seq, target_allele="", aligner=None):
        """ Align the seq on msa (Experimental)"""
        if not len(self.alleles):
            raise ValueError("MSA is empty")
        if not target_allele:
            target_allele = self.get_first()[0]
        if target_allele not in self.alleles:
            raise ValueError(f"{target_allele} not found")

        # setup aligner
        if not aligner:
            aligner = PairwiseAligner()
            aligner.alphabet = "ATCG-E"
            aligner.target_open_gap_score = -99999
            aligner.query_open_gap_score = -2

        # align
        target_seq = self.alleles[target_allele]
        result_seq = format(aligner.align(target_seq, seq)[0]).split("\n")[-2]
        assert len(result_seq) == len(target_seq)
        return result_seq

    # Functions for writing to bam file
    def _calculate_cigar(self, a, b):
        """ Compare two sequences and output cigar """
        # Compare two sequence
        c = []
        a = a.replace("E", "-")
        b = b.replace("E", "-")
        for i in range(len(a)):
            if a[i] == "-" and b[i] == "-":
                continue
            if a[i] == b[i]:
                c.append("M")
            elif a[i] == "-":
                c.append("I")
            elif b[i] == "-":
                c.append("D")
            else:
                c.append("X")

        # Aggregate the comparsion
        cigar = []
        for i in c:
            if not len(cigar):
                cigar.append([i, 1])
            elif cigar[-1][0] == i:
                cigar[-1][1] += 1
            else:
                cigar.append([i, 1])

        # Rename to cigar_tuple
        for i in cigar:
            if i[0] == "M":
                i[0] = 0
            elif i[0] == "I":
                i[0] = 1
            elif i[0] == "D":
                i[0] = 2
            elif i[0] == "X":
                i[0] = 8

        return cigar

    def save_bam(self, fname, ref_allele):
        """ Save the MSA to bam aligned on ref_allele """
        if not len(self.alleles):
            raise ValueError("MSA is empty")
        if ref_allele not in self.alleles:
            raise ValueError(f"{ref_allele} not found")
        if not fname:
            raise ValueError("filename is required")

        # setup reference and header
        if not ref_allele:
            ref_allele, ref_seq = self.get_first()
        else:
            ref_seq = self.alleles[ref_allele]
        header = {'HD': {'VN': '1.0'},
                  'SQ': [{'LN': len(ref_seq.replace("-", "").replace("E", "")),
                          'SN': ref_allele}]}

        # write bam file
        cigars = []
        with pysam.AlignmentFile(fname, "wb", header=header) as outf:
            for allele, seq in self.alleles.items():
                a = pysam.AlignedSegment()
                a.query_name = allele
                a.query_sequence = seq.replace("E", "").replace("-", "")
                cigars.append(self._calculate_cigar(ref_seq, seq))
                a.cigar = cigars[-1]

                # set to default
                a.mapping_quality = 60
                a.reference_id = 0
                a.reference_start = 0
                # a.template_length = 0
                # a.query_qualities = [30] * len(a.query_sequence)
                # a.flag = 0
                # a.next_reference_id = 0
                # a.next_reference_start = 0
                outf.write(a)
        pysam.sort("-o", fname, fname)
        pysam.index(fname)
        return cigars

    def _get_label(self):
        """ return assume label name for each chunk """
        if self.seq_type == "gen":
            assert len(self.blocks) % 2 == 1
            labels = []
            for i in range(len(self.blocks)):
                if i % 2:
                    labels.append(("exon", f"exon{i // 2 + 1}"))
                else:
                    labels.append(("intron", f"intron{i // 2}"))
            labels[0] = ("five_prime_UTR", "5UTR")
            labels[-1] = ("three_prime_UTR", "3UTR")
            return labels

        elif self.seq_type == "nuc":
            return [("exon", f"exon{i+1}") for i in range(len(self.blocks))]
        else:
            return [("gene_fragment", f"chunk{i+1}") for i in range(len(self.blocks))]

    def save_gff(self, fname, labels=[]):
        """ Output GFF3 format """
        # http://gmod.org/wiki/GFF3
        if not len(self.blocks):
            raise ValueError("MSA is empty")
        if not labels:
            labels = self._get_label()
        elif len(self.blocks) != len(labels):
            raise ValueError("Labels length is different from MSA chunks")

        block_pos = self.calculate_position()
        records = []

        # save allele info in each record
        for allele, seq in self.alleles.items():
            record = []
            pos = [len(seq[:i].replace("E", "").replace("-", "")) for i in block_pos]

            for i in range(len(self.blocks)):
                # ref source type start end . strand . tags
                record.append(
                    [allele, "pyHLAMSA", labels[i][0],
                     str(pos[i] + 1), str(pos[i + 1]), ".", "+", ".",
                     f"ID={labels[i][1]}_{allele}"]
                )
            records.extend(record)

        # save
        with open(fname, "w") as f:
            f.write("##gff-version 3\n")
            for record in records:
                f.write("\t".join(record) + "\n")

        return records


if __name__ == "__main__":
    # Basic operation: read, add
    hla = HLAmsa(["A", "B"], filetype=["gen", "nuc"], imgt_folder="alignments", version="3430")
    a = hla.genes["A"]
    a.add("A*consensus", a.get_consensus(include_gap=False))
    a.fill_imcomplete("A*consensus")

    # select and get consensus
    a_sub = a.select_allele(r"A\*.*:01:01:01$")
    a_sub.extend(a.select_allele(r"A\*consensus$"))
    print(a_sub)

    # output
    print(a_sub.select_exon().format_alignment_diff("A*consensus"))
    print(a_sub.to_MultipleSeqAlignment())
    a_sub.save_bam("tmp.bam", "A*consensus")
    SeqIO.write(a_sub.to_fasta(gap=False), "tmp.fa", "fasta")
    a_sub.save_gff("tmp.gff")

    # align seq on the consensus and print
    seq = list(a.get("A*01:01:87"))
    seq[55] = seq[65] = seq[78] = seq[79] = seq[80] = "T"
    seq = a_sub.align(''.join(seq), target_allele="A*consensus")
    a_sub.add("A*query", seq)
    print(a_sub[30:180].shrink().format_alignment_diff("A*consensus"))
