from __future__ import annotations
import re
import copy
import logging
from typing import List, Tuple, Dict, Tuple

from Bio.Align import MultipleSeqAlignment, PairwiseAligner
from Bio import AlignIO, SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import pysam


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

        blocks (list of dict):
            The dictionary contains three information

            * length: The length of each block
            * type: the name of block defined in Category:SO:SOFA. (Optional)
            * name: the name you want to show. (Optional)
    """

    def __init__(self, gene_name: str, seq_type="",
                 blocks=[]):
        self.gene_name = gene_name
        self.seq_type = seq_type
        self.alleles = {}
        self.blocks = copy.deepcopy(blocks)  # intron exon length
        self.logger = logging.getLogger(__name__)

    # Show the MSA attribute
    def __str__(self):
        return f"<{self.gene_name} {self.seq_type} "\
               f"alleles={len(self.alleles)} "\
               f"block={self.blocks}>"

    def get_length(self) -> int:
        """ Get the length of MSA """
        # 0 sequences is allow
        return sum(self.get_block_length())

    def get_block_length(self) -> int:
        """ Get the block's length of MSA """
        return [i['length'] for i in self.blocks]

    def get_sequence_num(self) -> int:
        """ Get the number of sequences in MSA """
        return len(self.alleles)

    def __len__(self) -> int:
        """ Get the number of sequences in MSA """
        return len(self.alleles)

    def size(self) -> Tuple[int, int]:
        """ Get the size (num_of_sequences, length_of_sequence) """
        return (len(self), self.get_length())

    def get_sequence_names(self) -> List[str]:
        """ Get the all the allele's sequence name in MSA """
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
                          self.blocks)
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

    def remove(self, name: Optional[str|List[str]]):
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
            new_msa.alleles = {allele: seq for allele, seq in self.alleles.items()
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
                i.count("-") + i.count("E")])
        return freqs

    def get_consensus(self, include_gap=False) -> str:
        """
        Generate the consensus sequence by choosing maximum frequency base

        Args:
          include_gap (bool):
              Allow consensus contains gap if gap is the maximum item.

              If include_gap=False and all the base on that position is gap (not shrinked before),
              it will warning and fill with A.
        """
        freqs = self.calculate_frequency()
        if not include_gap:
            if any(sum(f[:4]) == 0 for f in freqs):
                self.logger.warning("MSA contains gap, try .shrink() before .get_consensus()")
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
            new_msa.blocks[i]['length'] = sum(masks[gen_pos[i]:gen_pos[i+1]])
        assert sum(masks) == new_msa.get_length()

        # remove base in allele
        for allele, seq in self.alleles.items():
            new_msa.alleles[allele] = "".join([seq[i] for i in range(len(seq)) if masks[i]])

        return new_msa

    def get_variantion_base(self) -> List[str]:
        """
        Get the base positions where variation occurs

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
        """ Concat 2 MSA """
        if set(self.get_sequence_names()) != set(msa.get_sequence_names()):
            raise ValueError("Can not concat these two MSA because allele is different")
        new_msa = self.copy()
        new_msa.blocks.extend(copy.deepcopy(msa.blocks))
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
            new_msa.blocks = [{"length": len(new_msa._get_first()[1])}]
            return new_msa

        elif isinstance(index, tuple) or isinstance(index, list):
            new_msa = Genemsa(self.gene_name)
            new_msa.alleles = {allele: "".join([seq[i] for i in index])
                               for allele, seq in self.alleles.items()}
            new_msa.blocks = [{"length": len(new_msa._get_first()[1])}]
            return new_msa
        # Fail
        else:
            raise TypeError("Bad usage")

    # some helper functions
    def _get_block_position(self) -> List[int]:
        """ Calculate the start position of each block """
        pos = [0]
        for b in self.blocks:
            pos.append(pos[-1] + b['length'])
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
            if not (max(exon_index) <= len(self.blocks) // 2 and min(exon_index) > 0):
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
            new_seq = "".join([gen_seq[gen_pos[i]:gen_pos[i + 1]] for i in index])
            new_msa.alleles[allele] = new_seq

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
        """ Fill the `E` in exon-only sequences with ref_allele sequence """
        if ref_allele not in self.alleles:
            raise KeyError(f"{ref_allele} not found")

        ref_seq = self.alleles[ref_allele]
        for allele, seq in self.alleles.items():
            if "E" in seq:
                self.alleles[allele] = "".join([seq[i] if seq[i] != "E" else ref_seq[i]
                                                for i in range(len(seq))])
        return self

    def merge_exon(self, msa_nuc: Genemsa, replace_exon=True):
        """
        Merge nuc MSA into gen MSA

        It's allow that nuc MSA has new allele name than gen MSA,
        Genemsa will add the sequence in MSA, and the intron will fill by `E`

        If the exon part of gen MSA is differnet (e.g. less gapped) from nuc MSA,
        Genemsa will try to merge if it can

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
        if (self.gene_name != msa_nuc.gene_name or
                len(self.blocks) != len(msa_nuc.blocks) * 2 + 1):
            raise ValueError("Check object's name and block are correct")

        msas_gen = self.split()
        msas_nuc = msa_nuc.split()
        gen_names = set(self.get_sequence_names())
        nuc_names = set(msa_nuc.get_sequence_names())
        new_msa = Genemsa(self.gene_name, self.seq_type)
        new_msa.alleles = {name: "" for name in gen_names | nuc_names}

        for i_gen in range(len(self.blocks)):
            # intron -> fill with E
            if i_gen % 2 == 0:
                for name in nuc_names - gen_names:
                    msas_gen[i_gen].append(name, "E" * self.blocks[i_gen]['length'])
                new_msa += msas_gen[i_gen]
            # exon -> check before merge
            elif i_gen % 2 == 1:
                i_nuc = i_gen // 2
                # content change or length change
                if (msas_nuc[i_nuc].get_length() != msas_gen[i_gen].get_length()
                        or any(msas_nuc[i_nuc].get(name) != msas_gen[i_gen].get(name) for name in (nuc_names & gen_names))):
                    # check before merge
                    if len(gen_names - nuc_names):
                        raise ValueError("Some alleles doesn't exist in nuc MSA")
                    if any(msas_nuc[i_nuc].get(name).replace("-", "") !=
                           msas_gen[i_gen].get(name).replace("-", "")
                           for name in gen_names):
                        raise ValueError("Some exon sequences in gen MSA is not same as in nuc MSA")
                new_msa += msas_nuc[i_nuc]
        return new_msa

    # Format function
    def format_alignment_diff(self, ref_allele="", position_header=True) -> str:
        """
        Print the sequences of all alleles diff from `ref_allele` sequence.

        The format is similiar to IMGT alignment format.

        Returns:
          str: A formatted string

        Examples:
          >>> a = msa.select_allele(r"A\\*.*:01:01:01$").select_exon([6,7])
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
        return new_msa.format_alignment(position_header=position_header)

    def format_alignment(self, wrap=100, position_header=True) -> str:
        """
        Print the MSA

        Args:
          wrap (int): The max-length to wrap sequences

        Returns:
          str: A formatted string

        Examples:
          >>> a = msa.select_allele(r"A\\*.*:01:01:01$").select_exon([6,7])
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
        block_pos = self._get_block_position()
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
            if position_header:
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

    def format_alignment_from_center(self, pos: int, left=5, right=5) -> str:
        """
        Print all alleles sequences from the center of specific position

        Args:
          pos (int): The position wanted at the center
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
        output_str = ""
        prev_pos = max(pos - left, 0)
        output_str += f" {'gDNA':<18} {' ' * (pos - prev_pos)}{pos}\n"
        output_str += f" {' '   * 18} {' ' * (pos - prev_pos)}|\n"
        new_msa = self[prev_pos:pos+right]
        output_str += new_msa.format_alignment_diff(position_header=False)
        return output_str

    def format_variantion_base(self) -> str:
        """
        A handy function to show all the variation between the alleles

        Returns:
          str: A formatted string
        """
        bases = self.get_variantion_base()
        output_str = f"Total variantion: {len(bases)}\n"

        # merge if two variant are too close
        merged_bases = []
        right = 5
        for b in bases:
            if len(merged_bases) and merged_bases[-1][1] + right * 2 >= b:
                merged_bases[-1][1] = b
            else:
                merged_bases.append([b, b])

        # format the string
        for b_left, b_right in merged_bases:
            output_str += self.format_alignment_from_center(b_left, right=b_right - b_left + right)
        return output_str

    # Save or transform type/filetype
    def to_MultipleSeqAlignment(self) -> MultipleSeqAlignment:
        """
        Transfer this object to MultipleSeqAlignment(biopython)

        This operation will lost the information of intron's, exon's position and labels
        """
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
        It will automatically generate the label according on `seq_type`. (Inplace)
        """
        if self.seq_type == "gen":
            assert len(self.blocks) % 2 == 1 and len(self.blocks) >= 3
            for i in range(len(self.blocks)):
                if i % 2:
                    self.blocks[i]['type'] = "exon"
                    self.blocks[i]['name'] = f"exon{i // 2 + 1}"
                else:
                    self.blocks[i]['type'] = "intron"
                    self.blocks[i]['name'] = f"intron{i // 2}"

            self.blocks[0]['type'] = "five_prime_UTR"
            self.blocks[0]['name'] = f"5UTR"
            self.blocks[-1]['type'] = "three_prime_UTR"
            self.blocks[-1]['name'] = f"3UTR"
        elif self.seq_type == "nuc":
            for i in range(len(self.blocks)):
                self.blocks[i]['type'] = "exon"
                self.blocks[i]['name'] = f"exon{i+1}"
        else:
            for i in range(len(self.blocks)):
                self.blocks[i]['type'] = "gene_fragment"
                self.blocks[i]['name'] = f"block{i+1}"

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
        if not all(i.get("type") for i in self.blocks):
            self._assume_label()

        # TODO: should I save strand == '-' in model?
        # init pos and strand
        records = []
        block_pos = self._get_block_position()

        # save allele info in each record
        for allele, seq in self.alleles.items():
            record = []
            pos = [len(seq[:i].replace("E", "").replace("-", "")) for i in block_pos]

            for i in range(len(self.blocks)):
                # ref source type start end . strand . tags
                record.append(
                    [allele, "pyHLAMSA", self.blocks[i]['type'],
                     str(pos[i] + 1), str(pos[i + 1]), ".", strand, ".",
                     f"ID={self.blocks[i]['name']}_{allele}"]
                )
            records.extend(record)

        # save
        with open(fname, "w") as f:
            f.write("##gff-version 3\n")
            for record in records:
                f.write("\t".join(record) + "\n")
        return records

    # Read file function
    @classmethod
    def from_MultipleSeqAlignment(cls, bio_msa: MultipleSeqAlignment) -> Genemsa:
        """
        Transfer MultipleSeqAlignment instance(biopython) to Genemsa

        See more details in [biopython](https://biopython.org/docs/1.75/api/Bio.Align.html#Bio.Align.MultipleSeqAlignment)
        """
        new_msa = Genemsa("")
        new_msa.blocks = [{'length': bio_msa.get_alignment_length()}]
        for seq in bio_msa:
            new_msa.alleles[seq.id] = str(seq.seq)
        return new_msa

    @classmethod
    def read_MSF_file(cls, file_msf: str) -> Genemsa:
        """ Read .msf file """
        return Genemsa.from_MultipleSeqAlignment(AlignIO.read(file_msf, "msf"))

    @classmethod
    def read_dat_block(cls, file_dat: str):
        """
        Read block information in .dat file

        Args:
          file_dat (str): File path to xx.dat

        Returns:
          dat_object(dict): {allele_name: [[block_name(str), start(int), end(int)],]
              , where block_name is intron1, exon1, ...
              , and the start, end is the position without gap
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

    # out model save and load
    # TODO: rewrite to read/load json meta information
    @classmethod
    def load_msa(cls, file_fasta, file_gff) -> Genemsa:
        """ load this object to fasta and gff """
        # read
        msa = cls.from_MultipleSeqAlignment(AlignIO.read(file_fasta, "fasta"))
        msa.blocks = []

        # the blocks and labels from reference
        for row in open(file_gff):
            if not row.startswith("pyHLAMSA*consensus\t"):
                continue
            row = row.split("\t")
            msa.blocks.append({
                'length': int(row[4]) - int(row[3]) + 1,
                'type': row[2],
                'name': row[-1][3:].split("_")[0],
            })

        # check
        assert len(msa.alleles["pyHLAMSA*consensus"]) == msa.get_length()
        msa.remove("pyHLAMSA*consensus")
        return msa

    def save_msa(self, file_fasta, file_gff):
        """ Save this object to fasta and gff """
        # use pyHLAMSA*consensus as reference
        self.append("pyHLAMSA*consensus", self.get_consensus(include_gap=False))
        self.save_fasta(file_fasta, gap=True)
        self.save_gff(file_gff)
