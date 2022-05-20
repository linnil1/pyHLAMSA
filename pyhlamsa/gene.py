"""
Genemsa: core object in pyhlamsa

It can do many msa operation along with index (column position) and
block(intron exon region).

Basic usage:
``` python
from pyhlamsa import Genemsa
```
"""
from __future__ import annotations
import re
import copy
import logging
import dataclasses
from typing import List, Dict, Tuple, Iterator, Any
from collections.abc import Iterable

from Bio.Align import MultipleSeqAlignment, PairwiseAligner
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from .utils import cigar
from .utils.text import msa_to_string


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

        alleles (dict of str,str): MSA data.

            Allele name as the key and the sequence string as the value.

            The sequence has basic bases, "A", "T", "C", "G", "-" for gap,
            "E" stands for error (Mostly because some sequence has exon part only,
            so I fill the intron with E.

        blocks (list of BlockInfo): list of block information
        index (list of IndexInfo): list of index(position) information
        reference (str): The reference allele of the msa (Optional)
    """
    def __init__(self, gene_name: str, blocks=[], index=[], reference=None):
        self.gene_name = gene_name
        self.alleles = {}  # type: dict[str, str]
        self.blocks = copy.deepcopy(blocks)  # intron exon length
        self.logger = logging.getLogger(__name__)
        self.index = copy.deepcopy(index)
        self.reference = reference

    # Show the MSA attribute
    def __str__(self):
        block_info = " ".join([f"{b.name}({b.length})" for b in self.blocks])
        return f"<{self.gene_name} "\
               f"alleles={len(self.alleles)} "\
               f"block={block_info}>"

    def __repr__(self):
        return self.__str__()

    def get_length(self) -> int:
        """ Get the length of MSA """
        # 0 sequences is allow
        return sum(self.get_block_length())

    def get_block_length(self) -> List[int]:
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
        """ Same as list_allele_names """
        return list(self.alleles.keys())

    def list_alleles(self) -> List[str]:
        """
        List all the allele's sequence name in MSA

        Example:
           >>> a_gen.list_alleles()[:3]
           ['A*01:01:01:01', 'A*01:01:01:02N', 'A*01:01:01:03']
        """
        return list(self.alleles.keys())

    def get(self, allele: str) -> str:
        """ Get the sequence by allele name """
        return self.alleles[allele]

    def copy(self, copy_allele=True) -> Genemsa:
        """
        Clone a new MSA.

        Args:
            copy_allele (bool): copy the sequences
        """
        new_msa = Genemsa(self.gene_name,
                          blocks=self.blocks, index=self.index,
                          reference=self.reference)
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

    def remove(self, name: str | Iterable[str]):
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

    def select_allele(self, query: str | List[str]) -> Genemsa:
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

    def get_cigar(self, target_allele: str, ref_allele="") -> List[Tuple[str, int]]:
        """
        Get the cigar string of target_allele from ref_allele

        If ref_allele not set,
        it will automatically find the reference by get_reference

        Return:
          cigar(List[op_str, num]): The list of operator and number of bases
        Exmaple:
          cigar = [(M, 1), (X, 1), (D, 2), (M, 2)]
        """
        if target_allele not in self.alleles:
            raise KeyError(f"{target_allele} not found")
        if not ref_allele:
            ref_allele = self.get_reference()[0]
        if ref_allele not in self.alleles:
            raise KeyError(f"{ref_allele} not found")

        return cigar.calculate_cigar(self.alleles[ref_allele],
                                     self.alleles[target_allele])

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
            target_allele = self.get_reference()[0]
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
                Each items contains the number of ATCG and gap.
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

    def get_variantion_base(self) -> List[int]:
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
            new_msa = Genemsa(self.gene_name, reference=self.reference)
            new_msa.alleles = {allele: seq[index]
                               for allele, seq in self.alleles.items()}
            new_msa.blocks = [BlockInfo(length=len(new_msa.get_reference()[1]))]
            new_msa.index = copy.deepcopy(self.index[index])
            return new_msa
        elif isinstance(index, tuple) or isinstance(index, list):
            new_msa = Genemsa(self.gene_name, reference=self.reference)
            new_msa.alleles = {allele: "".join([seq[i] for i in index])
                               for allele, seq in self.alleles.items()}
            new_msa.blocks = [BlockInfo(length=len(new_msa.get_reference()[1]))]
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

    def set_reference(self, allele):
        """ Set the reference in msa (Inplace) """
        if allele not in self.alleles:
            raise IndexError(f"Cannot find {allele} in msa")
        self.reference = allele
        return self

    def get_reference(self) -> Tuple[str, str]:
        """
        Get the reference in MSA, if not existed, output the first one

        Returns:
          (allele_name, allele_seq)
        """
        if self.reference in self.alleles:
            return (self.reference, self.alleles[self.reference])
        if not len(self.alleles):
            raise ValueError("MSA is empty")
        return next(iter(self.alleles.items()))

    # Block-wise operation
    def select_exon(self, exon_index=[]) -> Genemsa:
        """
        Extract the exon by index.

        Args:
          exon_index (list of int or list of str): Index start from 1. i.e.

            * 1 for exon1
            * 2 for exon2

            Leave empty if you want all the exons

            If the exon_index contains list of string,
            it will select by name

        Example:
          >>> a_gen.select_exon([2]))
          <A gen alleles=7350 block=exon2(351)>

          >>> a_gen.select_exon([2, 3]))
          <A nuc alleles=7350 block=exon2(351) exon3(436)>

          >>> a_gen.select_exon(["exon2", "exon3"]))
          <A nuc alleles=7350 block=exon2(351) exon3(436)>
        """
        possible_exon_index = [i for i in range(len(self.blocks))
                               if self.blocks[i].type == "exon"]
        # If not specific the index, extract all exons
        if not exon_index:
            exon_index = possible_exon_index
        else:
            exon_indexs = []
            for i in exon_index:
                if type(i) is int:
                    exon_indexs.append(i * 2 - 1)
                else:
                    exon_indexs.append([b.name for b in self.blocks].index(i))
            exon_index = exon_indexs

        # check
        for ind in exon_index:
            if ind not in possible_exon_index:
                raise IndexError(f"You select the block is not exon: {ind}")
        new_msa = self.select_block(exon_index)
        return new_msa

    def select_block(self, index) -> Genemsa:
        """
        Extract blocks by index

        Args:
          index (list of int): Leave empty if you want all the blocks.

            Index start from 0.  e.g.

            * 0 for 5-UTR
            * 1 for exon1
            * 2 for intron1
            * 3 for exon2
            * 4 for 3-UTR(for two exons gene)
            * -1 for last block

            or you can use list of string,
            it will select by block name

        Example:
          >>> a_gen.select_block([-1])
          <A  alleles=4101 block=3UTR(302)>

          >>> a_gen.select_block([2, 3])
          <A  alleles=4101 block=intron1(130) exon2(335)>

          >>> a_gen.select_block(["5UTR", "exon3"])
          <A  alleles=4101 block=5UTR(301) exon3(413)>
        """
        block_name_mapping = {b.name: i for i, b in enumerate(self.blocks)}

        # replace -1 (3-UTR)
        for i in range(len(index)):
            if type(index[i]) is str:
                index[i] = block_name_mapping[index[i]]
            elif index[i] == -1:
                index[i] = len(self.blocks) - 1

        # check index boundary
        if not (max(index) < len(self.blocks) and min(index) >= -1):
            raise IndexError("Check block index is correct")

        # new a msa object
        new_msa = Genemsa(self.gene_name, reference=self.reference)
        for i in index:
            new_msa.blocks.append(copy.deepcopy(self.blocks[i]))

        # extract the sequences inside block region
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
        # A mapping from gen name to nuc index
        nuc_name_index = {b.name: i for i, b in enumerate(msa_nuc.blocks)
                          if b.type == "exon"}

        # check it's one-to-one mapping
        exon_set = set([b.name for b in self.blocks if b.type == "exon"])
        if set(nuc_name_index.keys()) != exon_set:
            raise ValueError(f"Cannot match blocks: "
                             f"gen={exon_set} nuc={nuc_name_index.keys()}")

        # create new msa and make sure the order of alleles
        new_msa = Genemsa(self.gene_name, reference=self.reference)
        new_msa.alleles = {name: "" for name in self.get_sequence_names()}
        new_msa.alleles.update({name: "" for name in msa_nuc.get_sequence_names()})

        # allele names
        gen_names = set(self.get_sequence_names())
        nuc_names = set(msa_nuc.get_sequence_names())
        exclude_name = set()  # type: set[str]

        # block-wise
        msas_gen = self.split()
        msas_nuc = msa_nuc.split()
        for i_gen in range(len(self.blocks)):
            # intron -> fill with E
            if self.blocks[i_gen].name not in nuc_name_index:
                for name in nuc_names - gen_names:
                    msas_gen[i_gen].append(name,
                                           "E" * self.blocks[i_gen].length)
                new_msa += msas_gen[i_gen].remove(list(exclude_name))
            # exon -> check before merge
            else:
                i_nuc = nuc_name_index[self.blocks[i_gen].name]
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
                    diff_names = list(diff_name)
                    if diff_names:
                        self.logger.warning(
                            f"Some exon sequences in gen MSA "
                            f"is not same as in nuc MSA "
                            f"{self.blocks[i_gen].name}: {diff_names}")
                        new_msa.remove(diff_names)
                        exclude_name.update(diff_names)
                new_msa += msas_nuc[i_nuc].remove(list(exclude_name))
        return new_msa.reset_index()

    # type transform
    def meta_to_json(self) -> Dict:
        """ Extract all meta information about this msa into json """
        meta = {
            'index': [dataclasses.asdict(i) for i in self.index],
            'blocks': [dataclasses.asdict(b) for b in self.blocks],
            'name': self.gene_name,
            'reference': self.reference,
        }
        return meta

    @classmethod
    def meta_from_json(cls, data: Dict[str, Any] = None) -> Genemsa:
        """ Import meta information from json """
        if data:
            return Genemsa(data['name'],
                           blocks=[BlockInfo(**b) for b in data['blocks']],
                           index=[IndexInfo(**i) for i in data['index']],
                           reference=data.get("reference"))
        else:
            return Genemsa("Unamed")

    def to_records(self: Genemsa, gap=True) -> list[SeqRecord]:
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

    def to_MultipleSeqAlignment(self) -> MultipleSeqAlignment:
        """ Transfer this object to MultipleSeqAlignment(biopython) """
        return MultipleSeqAlignment(self.to_records(gap=True),
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

    def assume_label(self, seq_type="gen") -> Genemsa:
        """
        It will automatically generate the block's label
        according on `seq_type`. (Inplace)

        seq_type:
          * gen: 5UTR-exon1-intron1-exon2-...-exon9-3UTR
          * nuc: exon1-exon2-...-exon9
          * other: block1-block2-block3-...
        """
        if seq_type == "gen":
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
        elif seq_type == "nuc":
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

    # Format function
    def format_alignment(self, wrap=100) -> str:
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
        if not self.index or not self.alleles:
            raise ValueError("MSA is empty")
        return msa_to_string(self, wrap=wrap)

    def format_alignment_diff(self, ref_allele="", show_position_set=None) -> str:
        """
        Print the sequences of all alleles diff from `ref_allele` sequence.

        The format is similiar to IMGT alignment format.

        * `-` indicate same as reference(first line)
        * `*` indicate deletion
        * `ATCG` indicate the SNV

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
        if not ref_allele:
            ref_allele = self.get_reference()[0]
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
        return msa_to_string(new_msa, show_position_set=show_position_set)

    def format_alignment_from_center(self, pos: int | Iterable[int],
                                     left=5, right=5) -> str:
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
            want_pos = [pos]
        else:
            want_pos = list(pos)  # type: ignore
        if not want_pos:
            return ""

        show_position_set = set(self.index[i].pos for i in want_pos)
        msa = None
        for p in want_pos:
            if msa is None:
                msa = self[p - left: p + right]
            else:
                msa += self[p - left: p + right]
        assert msa
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
        if not bases:
            return output_str

        # merge if two variant are too close
        merged_bases = []  # type: list[list[int]]
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
        assert msa
        return output_str + msa.format_alignment_diff(show_position_set=show_position_set)
