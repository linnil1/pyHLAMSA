from __future__ import annotations
import logging
import re
from pprint import pprint
from typing import List, Tuple, Dict

from Bio.Align import MultipleSeqAlignment, PairwiseAligner
from Bio import AlignIO, SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import pysam
import copy


class Genemsa:
    """
    An useful MSA interface

    Attributes:
        gene_name (str): The name of the gene
        seq_type (str): The type of the MSA.

            It can be `gen` and `nuc` or empty string if cannot determined

        alleles (dict of str,str): MSA data.

            Allele name as the key and the sequence string as the value

        blocks (list of int): The length of each chunk
        labels (list of str,str): The label of each chunk.

            Each tuple has two items: `label_type` and `label_name`.

            * `label_type` is the name of chunk defined in Category:SO:SOFA.
            * `label_name` is the name you want to show.
    """

    def __init__(self, gene_name: str, seq_type="",
                 blocks=[], labels=[]):
        self.gene_name = gene_name
        self.seq_type = seq_type
        self.alleles = {}
        self.blocks = copy.deepcopy(blocks)  # intron exon length
        self.labels = copy.deepcopy(labels)  # the label of the blocks
        self.logger = logging.getLogger(__name__)

    def __str__(self):
        return f"<{self.gene_name} {self.seq_type} "\
               f"alleles={len(self.alleles)} "\
               f"block={self.labels} length={self.blocks}>"

    # some helper functions
    def _calculate_position(self) -> List[int]:
        """ Calculate the start position of each chunk """
        pos = [0]
        for length in self.blocks:
            pos.append(pos[-1] + length)
        return pos

    def get_length(self) -> int:
        """ Get the length of MSA """
        return sum(self.blocks)

    def _get_first(self) -> Tuple[str, str]:
        """ Get the first record in MSA """
        return next(iter(self.alleles.items()))

    # reading functions
    def read_alignment_file(self, fname: str):
        """ Read MSA format file `fname` and save it in the instance """
        alleles = self.parse_alignment(fname)
        for allele, seq in alleles.items():
            self.alleles[allele] = seq.replace("|", "")
            # calculate length of exons and introns
            if not self.blocks:
                self.blocks = [len(seq) for seq in seq.split("|")]
        self.labels = self._assume_label()

    def parse_alignment(self, fname: str) -> Dict[str, str]:
        """
        Read MSA format file `fname`

        Returns:
          alleles (dict of str,str): The dictionary that map allele name to sequence
        """
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
            alleles[allele] = "".join(seq).replace(".", "-")

        for allele in rm_allele:
            self.logger.warning(f"Remove {allele} due to length in {fname}")
            del alleles[allele]
        return alleles

    # Position selection functions
    def select_exon(self, exon_index=[]) -> Genemsa:
        """
        Extract the exon by index.

        Args:
          exon_index (list of int): Index start from 1. i.e.

            * 1 for exon1
            * 2 for exon2

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

    def select_chunk(self, index=[]) -> Genemsa:
        """
        Extract chunks by index

        Args:
          index (list of int): Leave empty if you want all the chunks.

            Index start from 1.
            e.g.

            * 0 for 5-UTR
            * 1 for exon1
            * 2 for intron1
            * 3 for exon2
            * 4 for 3-UTR(for two exons gene)
            * -1 for 3-UTR(for all case)
        """
        if not index:
            return copy.deepcopy(self)

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
            new_msa.labels.append(self.labels[i])

        # extract
        gen_pos = self._calculate_position()
        for allele, gen_seq in self.alleles.items():
            new_seq = "".join([gen_seq[gen_pos[i]:gen_pos[i + 1]] for i in index])
            new_msa.alleles[allele] = new_seq

        return new_msa

    def __getitem__(self, index=[]) -> Genemsa:
        """
        Extract the region of the sequences by index

        Example:
          >>> msa = Genemsa("A", "gen")
          >>> msa.read_alignment_file("A_gen.txt")
          >>> # Inspect 50-100bp in the MSA
          >>> extract_msa = msa[50:100]
          >>> print(extract_msa)

          >>> # Inspect 2nd 3rd 5th bp in the MSA
          >>> extract_msa = msa[[1,2,4]]
          >>> print(extract_msa)
        """
        if not index:
            return self

        # Extract specific region in alignment
        if isinstance(index, slice) or isinstance(index, int):
            new_msa = Genemsa(self.gene_name)
            new_msa.alleles = {allele: seq[index]
                               for allele, seq in self.alleles.items()}
            new_msa.blocks = [len(new_msa._get_first()[1])]
            return new_msa

        elif isinstance(index, tuple) or isinstance(index, list):
            new_msa = Genemsa(self.gene_name)
            new_msa.blocks = [len(index)]
            new_msa.alleles = {allele: "".join([seq[i] for i in index])
                               for allele, seq in self.alleles.items()}
            return new_msa
        # Fail
        else:
            raise TypeError("Bad usage")

    # reverse
    def reverse_complement(self) -> Genemsa:
        """ Reverse the sequences """
        new_msa = Genemsa(self.gene_name, self.seq_type,
                          list(reversed(self.blocks)), self.labels)
        new_msa.alleles = {allele: str(Seq(seq).reverse_complement())
                           for allele, seq in self.alleles.items()}
        return new_msa

    # sort
    def sort(self):
        """ Sort the sequences """
        self.alleles = dict(sorted(self.alleles.items(), key=lambda i: i[1]))
        return self

    def sort_name(self):
        """ Sort the sequences by name """
        self.alleles = dict(sorted(self.alleles.items(), key=lambda i: i[0]))
        return self

    # Allele selection functions
    def get(self, ref_allele: str) -> str:
        """ Get the sequence by allele name """
        if ref_allele not in self.alleles:
            raise ValueError(f"{ref_allele} not found")
        return self.alleles[ref_allele]

    def select_allele(self, regex: str) -> Genemsa:
        """
        Select allele name by regex

        Args:
          regex (str): The regex pattern

        Examples:
          >>> # select allele name start with A*01:01
          >>> msa.select_allele("A\*01\:01.*")
        """
        new_msa = Genemsa(self.gene_name, self.seq_type,
                          self.blocks, self.labels)
        new_msa.alleles = {allele: seq for allele, seq in self.alleles.items()
                           if re.match(regex, allele)}
        return new_msa

    def select_complete(self) -> Genemsa:
        """ Select non exon-only sequences (No `E` in the sequence)"""
        new_msa = Genemsa(self.gene_name, self.seq_type,
                          self.blocks, self.labels)
        new_msa.alleles = {allele: seq for allele, seq in self.alleles.items()
                           if "E" not in seq}
        return new_msa

    def select_imcomplete(self) -> Genemsa:
        """ Select exon-only sequences (`E` exists in the sequence)"""
        new_msa = Genemsa(self.gene_name, self.seq_type,
                          self.blocks, self.labels)
        new_msa.alleles = {allele: seq for allele, seq in self.alleles.items()
                           if "E" in seq}
        return new_msa

    # Sequence adding functions
    def add(self, name: str, seq: str):
        """
        Add a sequence into MSA

        Make sure the sequence length is same as in MSA
        """
        if len(seq) != self.get_length():
            raise ValueError("Length not match to alignments")
        if not len(self.blocks):
            raise ValueError("MSA is empty")
        if name in self.alleles:
            raise ValueError(f"{name} already exist")

        self.alleles[name] = seq
        return self

    def extend(self, msa: Genemsa):
        """ Merge MSA into this instance """
        if self.blocks != msa.blocks:
            raise ValueError("Length is different")
        self.alleles.update(msa.alleles)
        return self

    # Functions deal with exon-only alleles
    def fill_imcomplete(self, ref_allele: str):
        """ Fill the `E` in exon-only sequences with ref_allele sequence """
        if ref_allele not in self.alleles:
            raise ValueError(f"{ref_allele} not found")

        ref_seq = self.alleles[ref_allele]
        for allele, seq in self.alleles.items():
            if "E" in seq:
                self.alleles[allele] = "".join([seq[i] if seq[i] != "E" else ref_seq[i]
                                                for i in range(len(seq))])
        return self

    def merge_exon(self, msa_nuc: Genemsa):
        """
        Merge nuc MSA into gen MSA

        The intron of nuc MSA will fill by `E`
        """
        # merge when allele in both nuc and gen
        if not (self.seq_type == "gen" and msa_nuc.seq_type == "nuc"):
            raise TypeError("Should merge nuc into gen")
        if (self.gene_name != msa_nuc.gene_name or
                len(self.blocks) != len(msa_nuc.blocks) * 2 + 1):
            raise ValueError("Check object's name and chunks are correct")

        nuc_pos = msa_nuc._calculate_position()
        gen_pos = self._calculate_position()
        new_alleles = {}
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
                            self.logger.warning(f"Remove {allele} nuc due to gen " \
                                           "is different from other gen after insert nuc")
                            break
                    new_seq += nuc_seq_part
                    blocks.append(len(nuc_seq_part))
                else:  # no break triggered
                    # check and add to new msa
                    if not ref_blocks:
                        ref_blocks = blocks
                    else:
                        assert ref_blocks == blocks
                    new_alleles[allele] = new_seq

            # Add exon only allele
            else:
                new_seq = ""
                for i in range(len(self.blocks)):
                    if i % 2 == 0:
                        new_seq += "E" * self.blocks[i]
                    else:
                        new_seq += nuc_seq[nuc_pos[i // 2]:nuc_pos[i // 2 + 1]]
                new_alleles[allele] = new_seq

        # Check
        non_insert_allele = set(self.alleles.keys()) - set(msa_nuc.alleles.keys())
        if len(non_insert_allele):
            if ref_blocks == self.blocks:
                for allele in non_insert_allele:
                    new_alleles[allele] = self.alleles[allele]
            else:
                self.logger.warning(f"Remove {non_insert_allele} due to " \
                               "gen is different from other gen after insert nuc")

        # overwrite
        self.blocks = ref_blocks
        self.alleles = new_alleles
        return self

    # Format function
    def format_alignment_diff(self, ref_allele="") -> str:
        """
        Print the sequences of all alleles diff from `ref_allele` sequence.

        The format is similiar to IMGT alignment format.

        Returns:
          str: A formatted string

        Examples:
          >>> a = msa.select_allele(r"A\*.*:01:01:01$").select_exon([6,7])
          >>> print(a.format_alignment_diff())
             gDNA               0
                                |
             A*01:01:01:01      ATAGAAAAGG AGGGAGTTAC ACTCAGGCTG CAA|GCAGTGA CAGTGCCCAG GGCTCTGATG TGTCTCTCAC AGCTTGTAAA G
             A*02:01:01:01      ---------- ------C--- T--------- ---|------- ---------- ---------- ---------- ---------- -
             A*03:01:01:01      ---------- ---------- ---------- ---|------- ---------- ---------- ----C----- ---------- -
             A*11:01:01:01      ---------- ---------- ---------- ---|------- ---------- ---------- ---------- ---------- -
             A*23:01:01:01      ---------- ------C--- T--------- ---|------- ---------- ---------- ---------- ---------- -
             A*25:01:01:01      ---------- ------C--- T--------- ---|------- ---------- ---------A ---------- ---------- -
             A*26:01:01:01      ---------- ------C--- T--------- ---|------- ---------- ---------A ---------- ---------- -
        """
        if not len(self.alleles):
            raise ValueError("MSA is empty")
        if not ref_allele:
            ref_allele = self._get_first()[0]
        if ref_allele not in self.alleles:
            raise ValueError(f"{ref_allele} not found")
        ref_seq = self.alleles[ref_allele]

        # use new msa object to save sequences
        new_msa = Genemsa(self.gene_name, self.seq_type,
                          self.blocks, self.labels)
        new_msa.alleles = {ref_allele: ref_seq}
        for allele, seq in self.alleles.items():
            if allele == ref_allele:
                continue
            new_seq = ""
            for i in range(len(seq)):
                if seq[i] == ref_seq[i]:
                    if ref_seq[i] == "-":
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

    def format_alignment(self, wrap=100) -> str:
        """
        Print the MSA

        Args:
          wrap (int): The max-length to wrap sequences

        Returns:
          str: A formatted string

        Examples:
          >>> a = msa.select_allele(r"A\*.*:01:01:01$").select_exon([6,7])
          >>> print(a.format_alignment())
             gDNA               0
                                |
             A*01:01:01:01      ATAGAAAAGG AGGGAGTTAC ACTCAGGCTG CAA|GCAGTGA CAGTGCCCAG GGCTCTGATG TGTCTCTCAC AGCTTGTAAA G
             A*02:01:01:01      ATAGAAAAGG AGGGAGCTAC TCTCAGGCTG CAA|GCAGTGA CAGTGCCCAG GGCTCTGATG TGTCTCTCAC AGCTTGTAAA G
             A*03:01:01:01      ATAGAAAAGG AGGGAGTTAC ACTCAGGCTG CAA|GCAGTGA CAGTGCCCAG GGCTCTGATG TGTCCCTCAC AGCTTGTAAA G
             A*11:01:01:01      ATAGAAAAGG AGGGAGTTAC ACTCAGGCTG CAA|GCAGTGA CAGTGCCCAG GGCTCTGATG TGTCTCTCAC AGCTTGTAAA G
             A*23:01:01:01      ATAGAAAAGG AGGGAGCTAC TCTCAGGCTG CAA|GCAGTGA CAGTGCCCAG GGCTCTGATG TGTCTCTCAC AGCTTGTAAA G
             A*25:01:01:01      ATAGAAAAGG AGGGAGCTAC TCTCAGGCTG CAA|GCAGTGA CAGTGCCCAG GGCTCTGATA TGTCTCTCAC AGCTTGTAAA G
             A*26:01:01:01      ATAGAAAAGG AGGGAGCTAC TCTCAGGCTG CAA|GCAGTGA CAGTGCCCAG GGCTCTGATA TGTCTCTCAC AGCTTGTAAA G
        """
        if not self.blocks or not self.alleles:
            raise ValueError("MSA is empty")
        output_str = ""
        bid = 0
        block_pos = self._calculate_position()
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
    def from_MultipleSeqAlignment(cls, bio_msa: MultipleSeqAlignment) -> Genemsa:
        """
        Transfer MultipleSeqAlignment instance(biopython) to Genemsa

        See more details in [biopython](https://biopython.org/docs/1.75/api/Bio.Align.html#Bio.Align.MultipleSeqAlignment)
        """
        new_msa = Genemsa("")
        new_msa.blocks = [bio_msa.get_alignment_length()]
        for seq in bio_msa:
            new_msa.alleles[seq.id] = str(seq.seq)
            if len(seq.seq) != new_msa.blocks[0]:
                raise ValueError("Length is different inside MultipleSeqAlignment")
        return new_msa

    def to_MultipleSeqAlignment(self) -> MultipleSeqAlignment:
        """ Transfer this object to MultipleSeqAlignment(biopython) """
        return MultipleSeqAlignment(self.to_fasta(gap=True))

    def to_fasta(self, gap=True) -> list[SeqRecord]:
        """
        Transfer MSA to list of SeqRecord

        Args:
            gap (bool): The sequence included gap or not
        """
        if gap:
            return [SeqRecord(Seq(seq.replace("E", "-")), id=allele, description="")
                    for allele, seq in self.alleles.items()]
        else:
            return [SeqRecord(Seq(seq.replace("E", "").replace("-", "")), id=allele, description="")
                    for allele, seq in self.alleles.items()]

    # Consensus function
    def calculate_frequency(self) -> List[List[int]]:
        """
        Calculate ATCG and gap frequency of each bp in MSA

        Returns:
            frequency (list of list of int):
                Each items contains the number of ATCG and space
        """
        freqs = []
        for i in zip(*self.alleles.values()):
            freqs.append([
                i.count("A"),
                i.count("T"),
                i.count("C"),
                i.count("G"),
                i.count("-")])
        return freqs

    def get_consensus(self, include_gap=False) -> str:
        """
        Generate the consensus sequence by choosing maximum frequency base

        Args:
          include_gap (bool): Allow consensus contains gap
        """
        freqs = self.calculate_frequency()
        if not include_gap:
            max_ind = [max(range(4), key=lambda i: f[i]) for f in freqs]
        else:
            max_ind = [max(range(5), key=lambda i: f[i]) for f in freqs]
        seq = ["ATCG-"[i] for i in max_ind]
        return "".join(seq)

    def shrink(self) -> Genemsa:
        """ Remove empty base """
        # index to delete
        freqs = self.calculate_frequency()
        masks = [f[4] != sum(f) for f in freqs]

        # recalcuate blocks
        new_msa = Genemsa(self.gene_name, self.seq_type,
                          [], self.labels)
        gen_pos = self._calculate_position()
        for i in range(len(self.blocks)):
            new_msa.blocks.append(sum(masks[gen_pos[i]:gen_pos[i+1]]))

        # remove bp in allele
        for allele, seq in self.alleles.items():
            new_msa.alleles[allele] = "".join([seq[i] for i in range(len(seq)) if masks[i]])
            assert len(new_msa.alleles[allele]) == new_msa.get_length()

        return new_msa

    # Align sequences on allele
    def align(self, seq: str, target_allele="", aligner=None) -> str:
        """
        Align the seq on msa (Experimental)

        Args:
            seq (str): The sequence you want to align on
            target_allele (str): I temporary align the sequences against target_allele
            aligner (PairwiseAligner): If set as None,
                the object will initizalize the parameters like EMBOSS
        Returns:
            result_sequence (str): The gap in sequence will set to '-'
        """
        if not len(self.alleles):
            raise ValueError("MSA is empty")
        if not target_allele:
            target_allele = self._get_first()[0]
        if target_allele not in self.alleles:
            raise ValueError(f"{target_allele} not found")

        # setup aligner
        if not aligner:
            aligner = PairwiseAligner()
            aligner.alphabet = "ATCG-E"
            aligner.target_open_gap_score = -99999
            aligner.query_open_gap_score = -10
            aligner.query_extend_gap_score = -.5
            aligner.match = 5
            aligner.mismatch = -4

        # align
        target_seq = self.alleles[target_allele]
        result_seq = format(aligner.align(target_seq, seq)[0]).split("\n")[-2]
        assert len(result_seq) == len(target_seq)
        return result_seq

    # Functions for writing to bam file
    def _calculate_cigar(self, a: str, b: str) -> List[Tuple[int, int]]:
        """
        Compare two sequences and output cigartuples

        The cigar_tuple is defined in https://pysam.readthedocs.io/en/latest/api.html#pysam.AlignedSegment.cigartuples
        """
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

    def save_fasta(self, fname: str, gap=True):
        """
        Save the MSA into fasta

        Args:
            gap (bool): The sequence included gap or not
        """
        SeqIO.write(self.to_fasta(gap=gap), fname, "fasta")

    def save_bam(self, fname: str, ref_allele="", save_ref=False):
        """
        Save the MSA into bam

        All alleles will seen as reads aligned on `ref_allele`

        Args:
          fname (str): The name of bamfile
          ref_allele (str): The reference allele.
              If the ref_allele is empty, the first allele will be reference.
          save_ref (bool): The reference allele will also save in the bam file
        """
        if not len(self.alleles):
            raise ValueError("MSA is empty")
        if not ref_allele:
            ref_allele = self._get_first()[0]
        if ref_allele not in self.alleles:
            raise ValueError(f"{ref_allele} not found")
        if not fname:
            raise ValueError("filename is required")

        # setup reference and header
        ref_seq = self.alleles[ref_allele]
        header = {'HD': {'VN': "1.0"},
                  'SQ': [{'LN': len(ref_seq.replace("-", "").replace("E", "")),
                          'SN': ref_allele}]}

        # write bam file
        cigars = []
        with pysam.AlignmentFile(fname, "wb", header=header) as outf:
            for allele, seq in self.alleles.items():
                # skip
                if not save_ref and allele == ref_allele:
                    continue

                # init bam record
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

    def _assume_label(self) -> List[str]:
        """
        Call this function when `self.labels` is empty.

        It will automatically generate the label according on `seq_type` and `blocks`.

        Returns:
          labels (list of str,str): The label of each chunks
        """
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

    def save_gff(self, fname: str, labels=[], strand="+"):
        """
        Save to GFF3 format

        Args:
          fname (str): The file name of gff3
          labels (list of str,str): Overwrite the labels
          strand (str): Must be "+" or "-".

              If the strand is "-", it will reverse the labels also.
        """
        # http://gmod.org/wiki/GFF3
        if not len(self.blocks):
            raise ValueError("MSA is empty")

        # labels
        if not labels:
            if not self.labels:
                self.labels = self._assume_label()
            labels = self.labels
            if strand == "-":
                labels = list(reversed(labels))
        elif len(self.blocks) != len(labels):
            raise ValueError("Labels length is different from MSA chunks")

        # init pos and strand
        records = []
        block_pos = self._calculate_position()

        # save allele info in each record
        for allele, seq in self.alleles.items():
            record = []
            pos = [len(seq[:i].replace("E", "").replace("-", "")) for i in block_pos]

            for i in range(len(self.blocks)):
                # ref source type start end . strand . tags
                record.append(
                    [allele, "pyHLAMSA", labels[i][0],
                     str(pos[i] + 1), str(pos[i + 1]), ".", strand, ".",
                     f"ID={labels[i][1]}_{allele}"]
                )
            records.extend(record)

        # save
        with open(fname, "w") as f:
            f.write("##gff-version 3\n")
            for record in records:
                f.write("\t".join(record) + "\n")
        return records

    # MSF related functions
    def read_MSF_file(self, file_msf):
        """ Read .msf file """
        new_msa = Genemsa.from_MultipleSeqAlignment(AlignIO.read(file_msf, "msf"))
        self.alleles = new_msa.alleles
        self.blocks = new_msa.blocks
        self.labels = new_msa.labels
        return self

    @classmethod
    def read_dat(cls, file_dat):
        """
        Read .dat file

        Args:
          file_dat (str): File path to xx.dat

        Returns:
          dat_object(dict): {allele_name: [[block_name(str), start(int), end(int)],]
              ,where block_name is intron1, exon1, ...
        """
        now_allele = ""
        read_next = False
        data = {}

        for line in open(file_dat):
            # read allele name
            if "allele=" in line:
                now_allele = line.split('"')[1]
                now_allele = now_allele.replace("HLA-", "")
                # no duplicated allele_name
                assert now_allele not in data
                data[now_allele] = []

            # read blocks
            elif re.match(r"FT\s+((UTR)|(exon)|(intron))", line):
                line = line.split()
                data[now_allele].append([line[1], *map(lambda a: int(a), line[2].split(".."))])
                if line[1] != "UTR":
                    read_next = True

            # read block name .e.g. exon1 exon2
            elif read_next:
                read_next = False
                data[now_allele][-1][0] += line.split('"')[1]

        # make sure the order start and end position in dat are correct
        for allele in data:
            data[allele] = list(sorted(data[allele], key=lambda a: a[1]))

        return data

    def merge_dat(self, dat) -> Genemsa:
        """
        Merge dat object

        Args:
          dat (object): Objected create from `Genemsa.read_dat`

        Returns:
            Genemsa
        """
        if not len(self.alleles):
            raise ValueError("MSA is empty")
        msf_length = len(self._get_first()[1])
        if self.seq_type != "gen" and self.seq_type != "nuc":
            raise TypeError("Check this object is gen or nuc")
        dat = copy.deepcopy(dat)

        # Check consistent between dat and msf
        alleles = self.alleles
        block_cord = {}
        new_alleles = {}
        for allele, seq in alleles.items():
            # Check allele name in the dat
            # and sequence is consistent to dat
            hla_name = allele
            assert hla_name in dat

            # extract exon region
            if self.seq_type == "nuc":
                dat_hla = []
                end = 0
                for i in [i for i in dat[hla_name] if "exon" in i[0]]:
                    gap = i[1] - end - 1
                    dat_hla.append([i[0], i[1] - gap, i[2] - gap])
                    end = dat_hla[-1][2]
                dat[hla_name] = dat_hla

            if len(seq.replace("-", "")) != dat[hla_name][-1][2]:
                self.logger.warning(f"Ignore {allele}, msf length is different from dat")
                continue

            # rename UTR
            # assume blocks are ordered
            if "UTR" == dat[hla_name][0][0]:
                dat[hla_name][0][0] = "5UTR"
            if "UTR" == dat[hla_name][-1][0]:
                dat[hla_name][-1][0] = "3UTR"

            # Coordination mapping from sequence to msf
            seq_cord = []
            assert len(seq) == msf_length
            for i in range(len(seq)):
                if seq[i] != "-":
                    assert seq[i] in "ATCG"
                    seq_cord.append(i)

            # Get the position of each block in msf
            old_cord = copy.deepcopy(block_cord)
            for block_name, start, end in dat[hla_name]:
                if block_name not in block_cord:
                    block_cord[block_name] = [99999, 0]

                block_cord[block_name][0] = min(block_cord[block_name][0],
                                                seq_cord[start - 1])
                block_cord[block_name][1] = max(block_cord[block_name][1],
                                                seq_cord[end - 1] + 1)

            # Check consistent
            # Assume first allele has correct coordination
            # , otherwise all allele will ignore
            block_cord_list = list(sorted(block_cord.values()))
            for i in range(1, len(block_cord_list)):
                if block_cord_list[i][0] < block_cord_list[i - 1][1]:
                    self.logger.warning(f"Ignore {allele}, msf inconsistent")
                    block_cord = old_cord
                    break
            else:
                new_alleles[allele] = seq

        # Reset position of each block
        # If perfect -> start_i == end_{i-1}
        # If not -> start_i = end_{i-1}
        block_cord_list = list(sorted(block_cord.items(), key=lambda i: i[1]))
        # TODO: is end == msf_length
        end = msf_length
        # end = block_cord_list[-1][1][1]
        for i in reversed(range(len(block_cord_list))):
            block_cord_list[i], end = ([block_cord_list[i][0],
                                        [block_cord_list[i][1][0], end]],
                                       block_cord_list[i][1][0])

        # Check the sequence length in msf is same as block coordination
        alleles = new_alleles
        new_alleles = {}
        block_cord = dict(block_cord_list)
        for allele, seq in alleles.items():
            allele_block = {i[0]: i[2] - i[1] + 1 for i in dat[allele]}

            for block_name in block_cord:
                len_block = len(seq[block_cord[block_name][0]:block_cord[block_name][1]].replace("-", ""))
                len_block_ref = allele_block.get(block_name, 0)
                # should always true, otherwise ignore in above
                assert len_block == len_block_ref
            else:
                new_alleles[allele] = seq

        # To new object
        new_msa = Genemsa(self.gene_name, self.seq_type)
        new_msa.alleles = new_alleles
        for name, (start, end) in block_cord_list:
            new_msa.blocks.append(end - start)
            if "3UTR" == name:
                new_msa.labels.append(("three_prime_UTR", "3UTR"))
            elif "5UTR" == name:
                new_msa.labels.append(("five_prime_UTR", "5UTR"))
            elif "exon" in name:
                new_msa.labels.append(("exon", name))
            elif "intron" in name:
                new_msa.labels.append(("intron", name))
        return new_msa

    @classmethod
    def load_msa(cls, file_fasta, file_gff) -> Genemsa:
        """ load this object to fasta and gff """
        # read
        msa = cls.from_MultipleSeqAlignment(AlignIO.read(file_fasta, "fasta"))
        msa.blocks = []
        msa.labels = []

        # the blocks and labels from reference
        for row in open(file_gff):
            if not row.startswith("pyHLAMSA*consensus\t"):
                continue
            row = row.split("\t")
            msa.labels.append([row[2], row[-1][3:].split("_")[0]])
            msa.blocks.append(int(row[4]) - int(row[3]) + 1)

        # check
        assert len(msa.alleles["pyHLAMSA*consensus"]) == sum(msa.blocks)
        del msa.alleles["pyHLAMSA*consensus"]
        return msa

    def save_msa(self, file_fasta, file_gff):
        """ Save this object to fasta and gff """
        # use pyHLAMSA*consensus as reference
        self.add("pyHLAMSA*consensus", self.get_consensus(include_gap=False))
        self.save_fasta(file_fasta, gap=True)
        self.save_gff(file_gff)
