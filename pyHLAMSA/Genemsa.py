from __future__ import annotations
import re
import copy
import json
import logging
import dataclasses
from typing import List, Tuple, Dict, Tuple

from Bio.Align import MultipleSeqAlignment, PairwiseAligner
from Bio import AlignIO, SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import pysam


@dataclasses.dataclass
class BlockInfo:
    """
    A class to save block information

    Attributes:
      length: The length of each block
      name: (Optional) The name of the block. e.g. intron1, exon2
      type: (Optional) The tag of the block defined in Category:SO:SOFA.
    """
    length: int
    name: str = ""
    type: str = ""


@dataclasses.dataclass
class IndexInfo:
    """
    A class to save index information

    Attributes:
      pos: The position index
      name: (Optional) The belonged block name for the position
      type: (Optional) The belonged block tag  for the position
    """
    pos: int
    name: str = ""
    type: str = ""


class Genemsa:
    """
    An useful MSA interface

    Attributes:
        gene_name (str): The name of the gene
        seq_type (str): The type of the MSA.

            It can be `gen` and `nuc` or empty string if cannot determined

        alleles (dict of str,str): MSA data.

            Allele name as the key and the sequence string as the value.

            The sequence has basic bases, "A", "T", "C", "G", "-" for gap,
            "E" stands for error (Mostly because some sequence has exon part only,
            so I fill the intron with E.

        blocks (list of BlockInfo): list of block information
        index (list of IndexInfo): list of index(position) information
    """

    def __init__(self, gene_name: str, seq_type="", blocks=[], index=[]):
        self.gene_name = gene_name
        self.seq_type = seq_type
        self.alleles = {}
        self.blocks = copy.deepcopy(blocks)  # intron exon length
        self.logger = logging.getLogger(__name__)
        self.index = copy.deepcopy(index)

    # Show the MSA attribute
    def __str__(self):
        block_info = ' '.join([f"{b.name}({b.length})" for b in self.blocks])
        return f"<{self.gene_name} {self.seq_type} "\
               f"alleles={len(self.alleles)} "\
               f"block={block_info}>"

    def get_length(self) -> int:
        """ Get the length of MSA """
        # 0 sequences is allow
        return sum(self.get_block_length())

    def get_block_length(self) -> int:
        """ Get the block's length of MSA """
        return [i.length for i in self.blocks]

    def __len__(self) -> int:
        """
        Get the number of sequences in MSA

        Example:
          >>> len(a_gen)
          4100
        """
        return len(self.alleles)

    def size(self) -> Tuple[int, int]:
        """ Get the size (num_of_sequences, length_of_sequence) """
        return (len(self), self.get_length())

    def get_sequence_names(self) -> List[str]:
        """
        Get the all the allele's sequence name in MSA

        Example:
           >>> a_gen.get_sequence_names()[:3]
           ['A*01:01:01:01', 'A*01:01:01:02N', 'A*01:01:01:03']
        """
        return list(self.alleles.keys())

    def get(self, ref_allele: str) -> str:
        """ Get the sequence by allele name """
        return self.alleles[ref_allele]

    def copy(self, copy_allele=True) -> Genemsa:
        """
        Clone a new MSA.

        Args:
            copy_allele (bool): copy the sequences
        """
        new_msa = Genemsa(self.gene_name, self.seq_type,
                          blocks=self.blocks, index=self.index)
        if copy_allele:
            new_msa.alleles = dict(self.alleles.items())
        return new_msa

    # allele-wise operation (row)
    def sort_name(self):
        """ Sort the sequences by name """
        self.alleles = dict(sorted(self.alleles.items(), key=lambda i: i[0]))
        return self

    def append(self, name: str, seq: str):
        """
        Add a sequence into MSA (inplace)

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

    def remove(self, name: Optional[str | List[str]]):
        """
        Remove a/some sequence from MSA (inplace)
        """
        if type(name) is str:
            del self.alleles[name]
            return self
        elif type(name) is list:
            for n in name:
                del self.alleles[n]
            return self
        # else
        raise NotImplementedError

    def select_allele(self, query: Optional[str, List[str]]) -> Genemsa:
        """
        Select allele name by regex or list of name

        Examples:
          >>> # select allele name start with A*01:01
          >>> msa.select_allele("A\\*01\\:01.*")
          >>> # select allele by list of string
          >>> msa.select_allele(["A*01:01", "A*02:03"])
        """
        new_msa = self.copy(copy_allele=False)
        if type(query) is str:
            new_msa.alleles = {allele: seq
                               for allele, seq in self.alleles.items()
                               if re.match(query, allele)}
        elif type(query) is list:
            new_msa.alleles = {name: self.alleles[name] for name in query}
        return new_msa

    def extend(self, msa: Genemsa):
        """ Merge MSA into this instance (inplace) """
        if self.get_block_length() != msa.get_block_length():
            raise ValueError("Length is different")
        leng = self.get_length()

        for name, seq in msa.alleles.items():
            if name in self.alleles:
                raise ValueError(f"{name} already exist")
            if len(seq) != leng:
                raise ValueError(f"Length is different, caused by {name}")
        self.alleles.update(msa.alleles)
        return self

    # sequence-level (row-level but consider bases)
    def sort(self):
        """ Sort the sequences """
        self.alleles = dict(sorted(self.alleles.items(), key=lambda i: i[1]))
        return self

    def reverse_complement(self) -> Genemsa:
        """ Reverse the sequences """
        new_msa = self.copy(copy_allele=False)
        new_msa.blocks = copy.deepcopy(list(reversed(self.blocks)))
        new_msa.index = copy.deepcopy(list(reversed(self.index)))
        new_msa.alleles = {allele: str(Seq(seq).reverse_complement())
                           for allele, seq in self.alleles.items()}
        return new_msa

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

    # base-wise operation (column)
    def calculate_frequency(self) -> List[List[int]]:
        """
        Calculate ATCG and gap frequency of each bp in MSA

        Returns:
            frequency (list of list of int):
                Each items contains the number of ATCG and gap. ("E" is count as gap)
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
          include_gap (bool):
              Allow consensus contains gap if gap is the maximum item.

              If include_gap=False and all the base on that position is gap
              (not shrinked before),
              it will warning and fill with A.
        """
        freqs = self.calculate_frequency()
        if not include_gap:
            if any(sum(f[:4]) == 0 for f in freqs):
                self.logger.warning(
                    "MSA contains gap, try .shrink() before .get_consensus()")
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
        new_msa = self.copy(copy_allele=False)

        # recalcuate blocks
        gen_pos = self._get_block_position()
        for i in range(len(self.blocks)):
            new_msa.blocks[i].length = sum(masks[gen_pos[i]:gen_pos[i + 1]])
        assert sum(masks) == new_msa.get_length()
        new_msa.index = [new_msa.index[i] for i in range(len(masks)) if masks[i]]

        # remove base in allele
        for allele, seq in self.alleles.items():
            new_msa.alleles[allele] = "".join(
                [seq[i] for i in range(len(seq)) if masks[i]])

        return new_msa

    def get_variantion_base(self) -> List[str]:
        """
        Get the base positions where variation occurs

        Example:
          ```
          msa:
            s0: AAT
            s1: AAC
            s2: CAC
          >>> msa.get_variantion_base()
          [0, 2]
          ```

        Returns:
          positions:
            Each integer represent the position of variation
        """
        freqs = self.calculate_frequency()
        num = len(self.alleles)
        base = []
        for i, freq in enumerate(freqs):
            if num not in freq:
                base.append(i)
        return base

    def __add__(self, msa: Genemsa) -> Genemsa:
        """
        Concat 2 MSA

        Example:
          >>> print(a_gen.select_exon([2]) + a_gen.select_exon([3]))
          <A nuc alleles=4100 block=exon2(335) exon3(413)>
        """
        names0 = set(self.get_sequence_names())
        names1 = set(msa.get_sequence_names())
        if names0 != names1:
            raise ValueError("Can not concat because some allele is miss: "
                             + str(names0.symmetric_difference(names1)))
        new_msa = self.copy()
        new_msa.blocks.extend(copy.deepcopy(msa.blocks))
        new_msa.index.extend(copy.deepcopy(msa.index))
        for name, seq in msa.alleles.items():
            new_msa.alleles[name] += seq
        return new_msa

    def __getitem__(self, index=[]) -> Genemsa:
        """
        Extract the region of the sequences by index (start from 0),
        but block information will not preserved

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
            new_msa.blocks = [BlockInfo(length=len(new_msa._get_first()[1]))]
            new_msa.index = copy.deepcopy(self.index[index])
            return new_msa
        elif isinstance(index, tuple) or isinstance(index, list):
            new_msa = Genemsa(self.gene_name)
            new_msa.alleles = {allele: "".join([seq[i] for i in index])
                               for allele, seq in self.alleles.items()}
            new_msa.blocks = [BlockInfo(length=len(new_msa._get_first()[1]))]
            new_msa.index = copy.deepcopy([self.index[i] for i in index])
            return new_msa
        # Fail
        else:
            raise TypeError("Bad usage")

    def reset_index(self) -> Genemsa:
        """
        Reset index:
        The old position information will be discard.

        Each position information will be counted from 0 and
        the label and name will copy from its block information
        """
        new_msa = self.copy()
        start = 0
        new_msa.index = []
        for b in self.blocks:
            for j in range(b.length):
                new_msa.index.append(IndexInfo(
                    pos=start,
                    type=b.type,
                    name=b.name,
                ))
                start += 1
        assert start == self.get_length()
        return new_msa

    # some helper functions
    def _get_block_position(self) -> List[int]:
        """ Calculate the start position of each block """
        pos = [0]
        for b in self.blocks:
            pos.append(pos[-1] + b.length)
        return pos

    def _get_first(self) -> Tuple[str, str]:
        """ Get the first record in MSA """
        return next(iter(self.alleles.items()))

    # Block-wise operation
    def select_exon(self, exon_index=[]) -> Genemsa:
        """
        Extract the exon by index.

        Args:
          exon_index (list of int): Index start from 1. i.e.

            * 1 for exon1
            * 2 for exon2

            Leave empty if you want all the exons
        """
        # TODO: use label to select
        if self.seq_type != "gen":
            raise TypeError("Check this object is gen or not")
        if len(self.blocks) % 2 != 1:
            raise IndexError("introns + exon should be odd")

        # If not specific the index, extract all exons
        if not exon_index:
            exon_index = range(1, len(self.blocks), 2)
        else:
            if not (max(exon_index) <= len(self.blocks) // 2
                    and min(exon_index) > 0):
                raise IndexError("Check block index is correct")
            exon_index = [i * 2 - 1 for i in exon_index]

        new_msa = self.select_block(exon_index)
        new_msa.seq_type = "nuc"
        # not need to check
        # for allele, seq in new_msa.alleles.items():
        #     assert "E" not in seq
        return new_msa

    def select_block(self, index=[]) -> Genemsa:
        """
        Extract blocks by index

        Args:
          index (list of int): Leave empty if you want all the blocks.

            Index start from 1.
            e.g.

            * 0 for 5-UTR
            * 1 for exon1
            * 2 for intron1
            * 3 for exon2
            * 4 for 3-UTR(for two exons gene)
            * -1 for 3-UTR(for all case)
        """
        if not index:  # all region
            return self.copy()
        if not (max(index) <= len(self.blocks) and min(index) >= -1):
            raise IndexError("Check block index is correct")

        # replace -1 (3-UTR)
        for i in range(len(index)):
            if index[i] == -1:
                index[i] = len(self.blocks) - 1

        # new a msa object
        new_msa = Genemsa(self.gene_name)
        for i in index:
            new_msa.blocks.append(copy.deepcopy(self.blocks[i]))

        # extract
        gen_pos = self._get_block_position()
        for allele, gen_seq in self.alleles.items():
            new_seq = "".join([gen_seq[gen_pos[i]:gen_pos[i + 1]]
                               for i in index])
            new_msa.alleles[allele] = new_seq

        # extract index
        for i in index:
            new_msa.index.extend(copy.deepcopy(self.index[gen_pos[i]:gen_pos[i + 1]]))
        return new_msa

    def split(self):
        """ Split the msa by blocks """
        return [self.select_block([i]) for i in range(len(self.blocks))]

    # Sequence-type operation(gen, nuc)
    def select_complete(self) -> Genemsa:
        """ Select non exon-only sequences (No `E` in the sequence)"""
        new_msa = self.copy(copy_allele=False)
        new_msa.alleles = {allele: seq for allele, seq in self.alleles.items()
                           if "E" not in seq}
        return new_msa

    def select_incomplete(self) -> Genemsa:
        """ Select exon-only sequences (`E` exists in the sequence)"""
        new_msa = self.copy(copy_allele=False)
        new_msa.alleles = {allele: seq for allele, seq in self.alleles.items()
                           if "E" in seq}
        return new_msa

    def fill_incomplete(self, ref_allele: str):
        """ Fill the `E` in exon-only sequences with ref_allele sequence (inplace) """
        if ref_allele not in self.alleles:
            raise KeyError(f"{ref_allele} not found")

        ref_seq = self.alleles[ref_allele]
        for allele, seq in self.alleles.items():
            if "E" in seq:
                self.alleles[allele] = "".join(
                    [seq[i] if seq[i] != "E" else ref_seq[i]
                     for i in range(len(seq))])
        return self

    def merge_exon(self, msa_nuc: Genemsa):
        """
        Merge nuc MSA into gen MSA

        It's allow that nuc MSA has new allele name than gen MSA,
        Genemsa will add the sequence in MSA, and the intron will fill by `E`

        If the exon part of gen MSA is differnet (e.g. less gapped) from nuc MSA,
        Genemsa will try to merge if it can

        Note that the index will be reset

        Example:
          ```
          # case1
          msa_gen:
            1: "AA|TT|CC",
            2: "AA|TC|CC",
          msa_nuc:
            3: "TC",
          After merge:
            1: "AA|TT|CC",
            2: "AA|TC|CC",
            3: "EE|TC|EE"
          ```

          ```
          # case2
          msa_gen:
            1: "AA|TT|CC",
            2: "AA|TC|CC",
          msa_nuc:
            1: "TT-",
            2: "T-C",
            4: "TTC",
          After merge:
            1: "AA|TT-|CC",
            2: "AA|T-C|CC",
            4: "EE|TTC|EE"
          ```
        """
        # TODO: rewrite: using label to merge
        if not (self.seq_type == "gen" and msa_nuc.seq_type == "nuc"):
            raise TypeError("Should merge nuc into gen")
        if (self.gene_name != msa_nuc.gene_name
                or len(self.blocks) != len(msa_nuc.blocks) * 2 + 1):
            raise ValueError("Check object's name and block are correct")

        msas_gen = self.split()
        msas_nuc = msa_nuc.split()
        gen_names = set(self.get_sequence_names())
        nuc_names = set(msa_nuc.get_sequence_names())
        new_msa = Genemsa(self.gene_name, self.seq_type)
        new_msa.alleles = {name: "" for name in gen_names | nuc_names}
        exclude_name = set()

        for i_gen in range(len(self.blocks)):
            # intron -> fill with E
            if i_gen % 2 == 0:
                for name in nuc_names - gen_names:
                    msas_gen[i_gen].append(name,
                                           "E" * self.blocks[i_gen].length)
                new_msa += msas_gen[i_gen].remove(list(exclude_name))
            # exon -> check before merge
            elif i_gen % 2 == 1:
                i_nuc = i_gen // 2
                # content change or length change
                if (msas_nuc[i_nuc].get_length() != msas_gen[i_gen].get_length()
                    or any(msas_nuc[i_nuc].get(name) != msas_gen[i_gen].get(name)
                           for name in (nuc_names & gen_names))):
                    # check before merge
                    if len(gen_names - nuc_names):
                        raise ValueError(
                            f"Some alleles doesn't exist in nuc MSA: "
                            f"{gen_names - nuc_names}")

                    diff_name = filter(lambda name:
                                       msas_nuc[i_nuc].get(name).replace("-", "")
                                       != msas_gen[i_gen].get(name).replace("-", ""),
                                       gen_names)
                    diff_name = list(diff_name)
                    if diff_name:
                        self.logger.warning(
                            f"Some exon sequences in gen MSA "
                            f"is not same as in nuc MSA "
                            f"{self.blocks[i_gen].name}: {diff_name}")
                        new_msa.remove(diff_name)
                        exclude_name.update(diff_name)
                new_msa += msas_nuc[i_nuc].remove(list(exclude_name))
        return new_msa.reset_index()

    # Format function
    def format_alignment_diff(self, ref_allele="", show_position_set=None) -> str:
        """
        Print the sequences of all alleles diff from `ref_allele` sequence.

        The format is similiar to IMGT alignment format.

        Check `.format_alignment()` for output detail.

        Returns:
          str: A formatted string

        Examples:
          >>> a = msa.select_allele(r"A\\*.*:01:01:01$").select_exon([6,7])
          >>> print(a.format_alignment_diff())
                              3166                                  3342
                                 |                                     |
             A*01:01:01:01       ATAGAAAAGG AGGGAGTTAC ACTCAGGCTG CAA| GCAGTGA CAGTGCC
             A*02:01:01:01       ---------- ------C--- T--------- ---| ------- -------
             A*03:01:01:01       ---------- ---------- ---------- ---| ------- -------
             A*11:01:01:01       ---------- ---------- ---------- ---| ------- -------
             A*23:01:01:01       ---------- ------C--- T--------- ---| ------- -------
             A*25:01:01:01       ---------- ------C--- T--------- ---| ------- -------
             A*26:01:01:01       ---------- ------C--- T--------- ---| ------- -------
             A*29:01:01:01       ---------- ------C--- T--------- ---| ------- -------
        """
        if not len(self.alleles):
            raise ValueError("MSA is empty")
        if not ref_allele:
            ref_allele = self._get_first()[0]
        if ref_allele not in self.alleles:
            raise ValueError(f"{ref_allele} not found")
        ref_seq = self.alleles[ref_allele]

        # use new msa object to save sequences
        new_msa = self.copy(copy_allele=False)
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
        return new_msa.format_alignment(show_position_set=show_position_set)

    def _apply_print_format(self, print_format_per_lines):
        """
        This function only controlling the layout not logic, used by `format_alignment`

        It may overlap the text if the print_format in `print_format_per_lines`
        doesn't given enough space.

        Args:
          print_format_per_lines (list of list of dict):
            The outer list is for each line and
            the inner list is a list of print_format
            Here is my print_format
            ```
            {
                # print the position
                index: True
                # type of printable word
                type: "base" | "allele_name" | "char"
                # the width
                length: 12  # for all
                # if type=base
                # pos = position of base
                pos: 12
                # if type=char
                # it will print the word
                word: ""  # for char
            }
            ```
        """
        output_str = ""
        for pfl in print_format_per_lines:
            # index line
            index_left_space = 0
            for pf in pfl:
                index_left_space += pf.get('length', 1)
                if pf.get('index'):
                    # note: 1-base
                    output_str += f"{self.index[pf['pos']].pos + 1:>{index_left_space}}"
                    index_left_space = 0
            output_str += "\n"

            # indicator line
            index_left_space = 0
            for pf in pfl:
                index_left_space += pf.get('length', 1)
                if pf.get('index'):
                    output_str += f"{'|':>{index_left_space}}"
                    index_left_space = 0
            output_str += "\n"

            # allele line
            for name, seq in self.alleles.items():
                for pf in pfl:
                    if pf['type'] == 'char':
                        output_str += pf['word']
                    elif pf['type'] == 'allele_name':
                        output_str += f"{name:<{pf['length']}}"
                    elif pf['type'] == 'base':
                        output_str += seq[pf['pos']]
                output_str += "\n"

            # separate each line
            output_str += "\n\n"
        return output_str

    def format_alignment(self, wrap=100, show_position_set=None) -> str:
        """
        Print the MSA

        Note: The index shown on output string is 1-base absolute position.
        If you want to print in relative position, use `.reset_index()` first.

        Args:
          wrap (int): The maximum base per line
          show_position_set (set): The list of position(0-base) you want to indicate

        Returns:
          str: A formatted string

        Examples:
          >>> a = msa.select_allele(r"A\\*.*:01:01:01$").select_exon([6,7])
          >>> print(a.format_alignment())
                              3166                                  3342
                                 |                                     |
             A*01:01:01:01       ATAGAAAAGG AGGGAGTTAC ACTCAGGCTG CAA| GCAGTGA CAGTGCCCAG
             A*02:01:01:01       ATAGAAAAGG AGGGAGCTAC TCTCAGGCTG CAA| GCAGTGA CAGTGCCCAG
             A*03:01:01:01       ATAGAAAAGG AGGGAGTTAC ACTCAGGCTG CAA| GCAGTGA CAGTGCCCAG
             A*11:01:01:01       ATAGAAAAGG AGGGAGTTAC ACTCAGGCTG CAA| GCAGTGA CAGTGCCCAG
             A*23:01:01:01       ATAGAAAAGG AGGGAGCTAC TCTCAGGCTG CAA| GCAGTGA CAGTGCCCAG
             A*25:01:01:01       ATAGAAAAGG AGGGAGCTAC TCTCAGGCTG CAA| GCAGTGA CAGTGCCCAG
        """
        if not self.blocks or not self.alleles:
            raise ValueError("MSA is empty")

        pos = 0
        last_index = -1
        seq_length = self.get_length()
        block_pos = self._get_block_position()[1:]  # remove first

        format_lines = []
        format_line = []
        while pos < seq_length:
            index = False  # if this base need index

            # init a line
            if not format_line:
                format_line.append({'type': 'char', 'word': " "})
                format_line.append({'type': 'allele_name', 'length': 18})
                format_line.append({'type': 'char', 'word': " "})
                line_length = 20
                last_index_pos = 0
                num_base_in_line = 0
                index = True

            # block indicator
            if pos == block_pos[0]:
                block_pos.pop(0)
                format_line.append({'type': 'char', 'word': "|"})
                index = True

            # not consecutive position: show the index
            if self.index[pos].pos != last_index + 1:
                format_line.append({'type': 'char', 'word': " "})
                line_length += 1
                index = True
            last_index = self.index[pos].pos

            # force to show
            if show_position_set is not None:
                index = self.index[pos].pos in show_position_set

            # space to avoid index overlapping
            # pos is 0-base
            while (index
                   and last_index_pos + len(str(self.index[pos].pos + 1)) > line_length):
                format_line.append({'type': 'char', 'word': " "})
                line_length += 1

            # base
            format_line.append({'type': 'base', 'pos': pos, 'index': index})
            pos += 1
            line_length += 1
            num_base_in_line += 1
            if index:
                last_index_pos = line_length

            # break the line
            if (num_base_in_line >= wrap
                    or (line_length > 120 and num_base_in_line % 10 == 0)):
                format_lines.append(format_line)
                format_line = []
                continue

            # split 10 bases in a line
            if num_base_in_line % 10 == 0:
                format_line.append({'type': 'char', 'word': " "})
                line_length += 1

        if format_line:
            format_lines.append(format_line)

        return self._apply_print_format(format_lines)

    def format_alignment_from_center(self, pos: int, left=5, right=5) -> str:
        """
        Print all alleles sequences from the center of specific position

        Check `.format_alignment()` for output format detail

        Args:
          pos (int or list of int): The (0-base relative) positions
          left (int): How many base shall print at the left of the center
          right (int): How many base shall print at the right of the center

        Returns:
          str: A formatted string

        Examples:
          ```
             gDNA                    3022
                                     |
             HLA-A*03:02        AGAGAAAAAT
             HLA-A*01:40        -----G----
             HLA-A*03:20        -----G----
          ```
        """
        if type(pos) is int:
            pos = [pos]

        show_position_set = set(self.index[i].pos for i in pos)
        msa = None
        for p in pos:
            if msa is None:
                msa = self[p - left: p + right]
            else:
                msa += self[p - left: p + right]
        return msa.format_alignment_diff(show_position_set=show_position_set)

    def format_variantion_base(self) -> str:
        """
        A handy function to show all the variation between the alleles

        Note: the `|` in the output string is NOT intron and exon boundary

        Check `.format_alignment()` for output format detail

        Returns:
          str: A formatted string
        Example:
            ```
            Total variantion: 71
                                    308    314          329          342    349
                                      |      |            |            |      |
             A*01:01:01:01       TGGCCGTCAT GGCGCC| CCCT CCTCCT| ACTC TCGGGGGCCC TGG
             A*02:01:01:01       ---------- ------| ---- -G----| ---- --------T- ---
             A*03:01:01:01       ---------- ------| ---- ------| ---- ---------- ---
             A*11:01:01:01       ---------- ------| ---- ------| ---- ---------- ---
             A*23:01:01:01       ---------- ------| ---- -G----| ---- ---------- ---
             A*25:01:01:01       ---------- ------| ---- -G----| ---- ---------- ---
             A*26:01:01:01       ---------- ------| ---- -G----| ---- ---------- ---
             A*29:01:01:01       ---------- ------| ---- ------| ---- -T-------- ---
             A*30:01:01:01       ---------- ------| ---- ------| ---- ---------- ---
             A*32:01:01:01       ---------- ------| ---- ------| ---- -T-------- ---
             A*33:01:01:01       ---------- ------| ---- ------| ---- -T-------- ---
             A*34:01:01:01       -----A---- ------| ---- -G----| ---- ---------- ---
             A*36:01:01:01       ---------- ------| ---- ------| ---- ---------- ---
             A*66:01:01:01       ---------- ------| ---- -G----| ---- ---------- ---
             A*68:01:01:01       ---------- ------| ---- -G----| ---- ---------- ---
             A*69:01:01:01       ---------- ------| ---- -G----| ---- ---------- ---
             A*74:01:01:01       ---------- ------| ---- ------| ---- -T-------- ---
             A*80:01:01:01       ---------- -C----| ---- ------| ---- ---------- ---
            ```
        """
        bases = self.get_variantion_base()
        output_str = f"#Total variantion: {len(bases)}\n"

        # merge if two variant are too close
        merged_bases = []
        right = 5
        for b in bases:
            if len(merged_bases) and merged_bases[-1][1] + right * 2 >= b:
                merged_bases[-1][1] = b
            else:
                merged_bases.append([b, b])

        # format the string
        # Using the property of index
        show_position_set = set(self.index[i].pos for i in bases)
        msa = None
        for b_left, b_right in merged_bases:
            if msa is None:
                msa = self[b_left - 5: b_right + 5]
            else:
                msa += self[b_left - 5: b_right + 5]
        return msa.format_alignment_diff(show_position_set=show_position_set)

    # Functions for writing to bam file
    def _calculate_cigar(self, a: str, b: str) -> List[Tuple[int, int]]:
        """
        Compare two sequences and output cigartuples

        The cigar_tuple is defined in
        https://pysam.readthedocs.io/en/latest/api.html#pysam.AlignedSegment.cigartuples
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

    # Save function
    def to_fasta(self, gap=True) -> list[SeqRecord]:
        """
        Transfer MSA to list of SeqRecord

        Args:
            gap (bool): The sequence included gap or not
        """
        if gap:
            return [SeqRecord(Seq(seq.replace("E", "-")),
                              id=allele, description="")
                    for allele, seq in self.alleles.items()]
        else:
            return [SeqRecord(Seq(seq.replace("E", "").replace("-", "")),
                              id=allele, description="")
                    for allele, seq in self.alleles.items()]

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
        It will automatically generate the label according on `seq_type`. (Inplace)
        """
        if self.seq_type == "gen":
            assert len(self.blocks) % 2 == 1 and len(self.blocks) >= 3
            for i in range(len(self.blocks)):
                if i % 2:
                    self.blocks[i].type = "exon"
                    self.blocks[i].name = f"exon{i // 2 + 1}"
                else:
                    self.blocks[i].type = "intron"
                    self.blocks[i].name = f"intron{i // 2}"

            self.blocks[0].type = "five_prime_UTR"
            self.blocks[0].name = f"5UTR"
            self.blocks[-1].type = "three_prime_UTR"
            self.blocks[-1].name = f"3UTR"
        elif self.seq_type == "nuc":
            for i in range(len(self.blocks)):
                self.blocks[i].type = "exon"
                self.blocks[i].name = f"exon{i+1}"
        else:
            for i in range(len(self.blocks)):
                self.blocks[i].type = "gene_fragment"
                self.blocks[i].name = f"block{i+1}"

        # inplace to reset the index
        self.index = self.reset_index().index
        return self

    def save_gff(self, fname: str, strand="+"):
        """
        Save to GFF3 format

        Args:
          fname (str): The file name of gff3
          strand (str): Must be "+" or "-".

              If the strand is "-", it will add `strand` in GFF file,
              if you want to reverse the sequence, please use `reverse_complement` first.
        """
        # http://gmod.org/wiki/GFF3
        if not len(self.blocks):
            raise ValueError("MSA is empty")

        # labels
        if not all(b.type for b in self.blocks):
            self._assume_label()

        # TODO: should I save strand == '-' in model?
        # init pos and strand
        records = []
        block_pos = self._get_block_position()

        # save allele info in each record
        for allele, seq in self.alleles.items():
            record = []
            pos = [len(seq[:i].replace("E", "").replace("-", ""))
                   for i in block_pos]

            for i in range(len(self.blocks)):
                # ref source type start end . strand . tags
                record.append(
                    [allele, "pyHLAMSA", self.blocks[i].type,
                     str(pos[i] + 1), str(pos[i + 1]), ".", strand, ".",
                     f"ID={self.blocks[i].name}_{allele}"]
                )
            records.extend(record)

        # save
        with open(fname, "w") as f:
            f.write("##gff-version 3\n")
            for record in records:
                f.write("\t".join(record) + "\n")
        return records

    # type transform
    def meta_to_json(self) -> Dict:
        """ Extract all meta information about this msa into json """
        meta = {
            'index': [dataclasses.asdict(i) for i in self.index],
            'blocks': [dataclasses.asdict(b) for b in self.blocks],
            'seq_type': self.seq_type,
            'name': self.gene_name,
        }
        return meta

    @classmethod
    def meta_from_json(cls, data=None) -> Genemsa:
        """ Import meta information from json """
        if data:
            return Genemsa(data['name'], data['seq_type'],
                           blocks=[BlockInfo(**b) for b in data['blocks']],
                           index=[IndexInfo(**i) for i in data['index']])
        else:
            return Genemsa("")

    def to_MultipleSeqAlignment(self) -> MultipleSeqAlignment:
        """ Transfer this object to MultipleSeqAlignment(biopython) """
        return MultipleSeqAlignment(self.to_fasta(gap=True),
                                    annotations=self.meta_to_json())

    @classmethod
    def from_MultipleSeqAlignment(cls, bio_msa: MultipleSeqAlignment) -> Genemsa:
        """
        Transfer MultipleSeqAlignment instance(biopython) to Genemsa

        See more details in [biopython](
        https://biopython.org/docs/1.75/api/Bio.Align.html#Bio.Align.MultipleSeqAlignment)
        """
        new_msa = cls.meta_from_json(bio_msa.annotations)
        if not new_msa.blocks:
            new_msa.blocks = [BlockInfo(length=bio_msa.get_alignment_length())]
        if not new_msa.index:
            new_msa = new_msa.reset_index()

        for seq in bio_msa:
            new_msa.alleles[seq.id] = str(seq.seq)
        assert new_msa.get_length() == bio_msa.get_alignment_length()
        return new_msa

    # out model save and load
    @classmethod
    def load_msa(cls, file_fasta, file_json) -> Genemsa:
        """
        load this object to fasta and gff

        Example:
          >>> from pyHLAMSA import Genemsa
          >>> a_gen.save_msa("a_gen.fa", "a_gen.json")
          >>> a_gen = Genemsa.load_msa("a_gen.fa", "a_gen.json")
        """
        # read
        msa = AlignIO.read(file_fasta, "fasta")
        with open(file_json) as f:
            msa.annotations = json.load(f)
        # main
        new_msa = cls.from_MultipleSeqAlignment(msa)
        assert len(new_msa._get_first()[1]) == new_msa.get_length()
        return new_msa

    def save_msa(self, file_fasta, file_json):
        """ Save this object to fasta and gff """
        self.save_fasta(file_fasta, gap=True)
        with open(file_json, "w") as f:
            json.dump(self.meta_to_json(), f)
