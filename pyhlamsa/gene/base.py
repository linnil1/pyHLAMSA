"""
Define the basic class of Genemsa,
including GenemsaBase, BlockInfo, IndexInfo
"""

import re
import copy
import logging
import dataclasses
from typing import TypeVar, Optional, Union, Any
from collections.abc import Iterable, KeysView, Iterator

from Bio.Seq import Seq

from ..utils import cigar


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


GenemsaType = TypeVar("GenemsaType", bound="GenemsaBase")
BlockInput = Union[int, str, BlockInfo]


class GenemsaBase:
    """
    The MSA has additional block-level infomations,
    which is useful for integrating dna and rna data.

    ```
    Block-level:
    name           intron1  exon1  intron2
    region         <-------><-----><------------>

    Index-level:
    (column-level) 123456789012345678901234567890

    Sequence-level:
    allele_name1   ATTTTTCTTTTTTGTTTTTTATTTTTTCTT
    allele_name2   ATTTTTCTTTTTTGTTTTTTATTTTTTCTT

    The sequence has basic bases, "A", "T", "C", "G" and "-" for gap,
    "E" stands for error.
    (Mostly because some sequence has exon part only, so I fill the intron with E.)
    ```
    """

    def __init__(
        self,
        gene_name: str = "Unamed",
        blocks: Iterable[BlockInfo] = [],
        index: Iterable[IndexInfo] = [],
        reference: str = "",
    ):
        """
        Attributes:
            gene_name: The name of the gene

            alleles MSA data.

                Allele name as the key and the sequence string as the value.


            blocks: list of block information
            index: list of index(position) information
            reference: allele of the msa (Optional)
        """
        self.gene_name = gene_name
        self.alleles: dict[str, str] = {}
        self.blocks = copy.deepcopy(list(blocks or []))  # intron exon length
        self.index = copy.deepcopy(list(index or []))
        self.logger = logging.getLogger(__name__)
        self.reference = reference

    def copy(self: GenemsaType, copy_allele: bool = True) -> GenemsaType:
        """
        Clone a new MSA.

        Args:
          copy_allele: Copy the sequences as well
        """
        # Child's Type
        Genemsa = type(self)
        new_msa = Genemsa(
            self.gene_name,
            blocks=self.blocks,
            index=self.index,
            reference=self.reference,
        )
        if copy_allele:
            new_msa.alleles = dict(self.items())
        return new_msa

    def __repr__(self) -> str:
        """Show name, number of alleles and block infos in this MSA"""
        block_info = " ".join([f"{b.name}({b.length})" for b in self.blocks])
        return f"<{self.gene_name} alleles={len(self)} block={block_info}>"

    # size
    def get_length(self) -> int:
        """Get the length of MSA"""
        # 0 sequences is allow
        if self:
            return len(self.get_reference()[1])
        else:
            # another way to calculate length
            return sum(i.length for i in self.blocks)

    def size(self) -> tuple[int, int]:
        """Get the size (num_of_allele, length_of_sequence)"""
        return (len(self), self.get_length())

    def list_alleles(self) -> KeysView[str]:
        """
        List all the allele's sequence name in MSA

        Example:
           >>> a_gen.list_alleles()[:3]
           ['A*01:01:01:01', 'A*01:01:01:02N', 'A*01:01:01:03']
        """
        return self.alleles.keys()

    def get_sequence_names(self) -> KeysView[str]:
        """Same as list_alleles"""
        return self.list_alleles()

    # dict-like function
    def __iter__(self) -> Iterator[str]:
        """Iter allele like iter(dict)"""
        return iter(self.alleles)

    def items(self) -> Iterable[tuple[str, str]]:
        """list allele name along with allele sequence like dict.items()"""
        return self.alleles.items()

    def get(self, allele: str) -> str:
        """Get the sequence by allele name"""
        return self.alleles[allele]

    def contains(self, allele: str) -> bool:
        """Implement `in` operator"""
        return allele in self.alleles

    def __len__(self) -> int:
        """
        Get the number of alleles in the MSA

        Example:
          >>> len(a_gen)
          4100
        """
        return len(self.alleles)

    def truth(self) -> bool:
        """
        If msa has 0 alleles return False
        else return Implement if(self) function
        """
        return bool(self.alleles)

    # reference-related function
    def set_reference(self: GenemsaType, allele: str) -> GenemsaType:
        """Set the reference in msa (Inplace)"""
        if allele not in self:
            self.logger.warning(f"Cannot find {allele} in MSA")
        self.reference = allele
        return self

    def get_reference(self) -> tuple[str, str]:
        """
        Get the reference in MSA, if not existed, output the first one

        Returns:
          (allele_name, allele_seq)
        """
        if self.reference in self:
            return (self.reference, self.get(self.reference))
        if not self:
            raise ValueError("MSA is empty")
        self.logger.warning(
            f"Reference {self.reference} not existed in MSA."
            " Using the first allele instead."
        )
        allele, seq = next(iter(self.items()))
        self.reference = allele  # set it
        return allele, seq

    def get_allele_or_error(self, allele: str = "") -> tuple[str, str]:
        """
        Get the sequence by allele name

        Args:
          allele (str): Allele name. If not provided, reference allele are used

        Returns:
          allele_name (str) and its sequence (str):
        """
        if not allele:
            allele = self.get_reference()[0]
        if allele not in self:
            raise ValueError(f"{allele} not found")
        return allele, self.get(allele)

    # allele function
    def sort_name(self: GenemsaType) -> GenemsaType:
        """Sort the allele by alelle name"""
        self.alleles = dict(sorted(self.items(), key=lambda i: i[0]))
        return self

    def append(self: GenemsaType, name: str, seq: str) -> GenemsaType:
        """
        Add a sequence into MSA (inplace)

        Make sure the sequence length is same as in MSA
        """
        if len(seq) != self.get_length():
            raise ValueError("Length not match to alignments")
        if not self.blocks:
            raise ValueError("MSA doesn't has block info")
        if name in self:
            raise ValueError(f"{name} already exist")

        self.alleles[name] = seq
        return self

    def extend(self: GenemsaType, msa: GenemsaType) -> GenemsaType:
        """Add MSA's alleles into this MSA (inplace)"""
        if [b.length for b in self.blocks] != [b.length for b in msa.blocks]:
            raise ValueError("Length is different")
        leng = self.get_length()

        for name, seq in msa.items():
            if name in self:
                raise ValueError(f"{name} already exist")
            if len(seq) != leng:
                raise ValueError(f"Length is different, caused by {name}")
        self.alleles.update(msa.alleles)
        return self

    def remove(self: GenemsaType, name: Union[str, Iterable[str]]) -> GenemsaType:
        """
        Remove a/some sequence from MSA (inplace)
        """
        if isinstance(name, str):
            del self.alleles[name]
        elif isinstance(name, (list, tuple)):
            for n in name:
                del self.alleles[n]
        return self
        # else
        raise NotImplementedError

    def select_allele(
        self: GenemsaType, query: Union[str, Iterable[str]]
    ) -> GenemsaType:
        """
        Select allele name by regex or list of name

        Examples:
          >>> # select allele name start with A*01:01
          >>> msa.select_allele(r"A\\*01:01.*")
          >>> # select allele by list of string
          >>> msa.select_allele(["A*01:01", "A*02:03"])
        """
        new_msa = self.copy(copy_allele=False)
        if isinstance(query, str):
            new_msa.alleles = {
                allele: seq for allele, seq in self.items() if re.match(query, allele)
            }
        elif isinstance(query, Iterable):
            new_msa.alleles = {name: self.alleles[name] for name in query}
        return new_msa

    # sequence functions
    def sort(self: GenemsaType) -> GenemsaType:
        """Sort the allele by their sequences"""
        self.alleles = dict(sorted(self.items(), key=lambda i: i[1]))
        return self

    def reverse_complement(self: GenemsaType) -> GenemsaType:
        """Reverse the sequences"""
        new_msa = self.copy(copy_allele=False)
        new_msa.blocks = copy.deepcopy(list(reversed(self.blocks)))
        new_msa.index = copy.deepcopy(list(reversed(self.index)))
        new_msa.alleles = {
            allele: str(Seq(seq).reverse_complement()) for allele, seq in self.items()
        }
        return new_msa

    def get_cigar(
        self, target_allele: str, ref_allele: str = ""
    ) -> list[tuple[str, int]]:
        """
        Get the cigar string of target_allele from ref_allele

        If ref_allele not set,
        it will automatically find the reference by get_reference

        Return:
          cigar(list[op_str, num]): The list of operator and number of bases
        Exmaple:
          `cigar = [(M, 1), (X, 1), (D, 2), (M, 2)]`
        """
        if target_allele not in self:
            raise KeyError(f"{target_allele} not found")
        ref_allele, _ = self.get_allele_or_error(ref_allele)
        return cigar.calculate_cigar(
            self.alleles[ref_allele], self.alleles[target_allele]
        )

    # column operation part
    def calculate_frequency(self) -> list[list[int]]:
        """
        Calculate ATCG and gap frequency of each bp in MSA

        Returns:
          frequency (list[list[int]]):
            Each items contains the number of ATCG and gap.
        """
        freqs = []
        for i in zip(*self.alleles.values()):
            freqs.append(
                [i.count("A"), i.count("T"), i.count("C"), i.count("G"), i.count("-")]
            )
        return freqs

    def get_consensus(self, include_gap: bool = False) -> str:
        """
        Generate the consensus sequence by choosing maximum frequency base

        Args:
          include_gap (bool):
            Allow consensus contains gap if gap is the maximum item.

            If include_gap=False and all the base on that position is gap
            (not shrinked before),
            it will warning and fill with A.

            `E` will be ignored.

        Example:
        ```
        a0                 CCATT-GGT--GTCGGGTTTCCAG
        a1                 CCACTGGGT--ATCGGGTTTCCAG
        c2                 CAATTGGGT--GTCGGGT---AAG
        consensus          CCATTGGGT--GTCGGGTTTCCAG
        consensus(no-gap)  CCATTGGGTAAGTCGGGTTTCCAG
        ```
        """
        freqs = self.calculate_frequency()
        if not include_gap:
            if any(sum(f[:4]) == 0 for f in freqs):
                self.logger.warning(
                    "MSA contains gap, try .shrink() before .get_consensus()"
                )
            max_ind = [max(range(4), key=lambda i: f[i]) for f in freqs]
        else:
            max_ind = [max(range(5), key=lambda i: f[i]) for f in freqs]
        return "".join(map(lambda i: "ATCG-"[i], max_ind))

    def get_variantion_base(self) -> list[int]:
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
        num = len(self)
        base = []
        for i, freq in enumerate(freqs):
            if num not in freq:
                base.append(i)
        return base

    # block functions
    def _get_block_index(self, block: BlockInput) -> int:
        """Find the index of the block"""
        if isinstance(block, str):
            for i, b in enumerate(self.blocks):
                if b.name == block:
                    return i
            raise ValueError(f"{block} not found")
        elif isinstance(block, BlockInfo):
            for i, b in enumerate(self.blocks):
                if b.name == block.name:
                    return i
            raise ValueError(f"{block} not found")
        elif isinstance(block, int):
            id = block
            if id < 0:
                id = len(self.blocks) + id
            if not (0 <= id < len(self.blocks)):
                raise IndexError(f"{block} is out of index")
            return id

        else:
            raise NotImplementedError(f"Type of {block} not work now")

    def get_block_interval(self, block: BlockInput) -> tuple[int, int]:
        """Calculate the start(included) and end index (excluded) of the block"""
        index = self._get_block_index(block)
        start = sum(self.blocks[i].length for i in range(index))
        return start, start + self.blocks[index].length

    def get_block_position(self, block: BlockInput) -> int:
        """Calculate the start position of the block"""
        index = self._get_block_index(block)
        return sum(self.blocks[i].length for i in range(index))

    def select_block(self: GenemsaType, index: Iterable[BlockInput]) -> GenemsaType:
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
        new_msa = self.copy(copy_allele=False)

        # choose block index by index
        new_block = []
        new_index = []
        all_pos = []
        index_ids = [self._get_block_index(i) for i in index]
        for i in index_ids:
            new_block.append(self.blocks[i])
            start, end = self.get_block_interval(i)
            all_pos.append((start, end))
            new_index.extend(self.index[start:end])
        new_msa.blocks = new_block
        new_msa.index = new_index

        # extract the sequences inside block region
        for allele, gen_seq in self.items():
            new_seq = "".join([gen_seq[start:end] for start, end in all_pos])
            new_msa.alleles[allele] = new_seq
        return new_msa.copy()

    def select_exon(
        self: GenemsaType, exon_index: Iterable[BlockInput] = []
    ) -> GenemsaType:
        """
        Extract the exon by index.

        Args:
          exon_index (list[str|int]): Index start from 1. i.e.

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
        exons = [b for b in self.blocks if b.type == "exon"]

        # If not specific the index, extract all exons
        exon_list: list[BlockInput] = []
        if not exon_index:
            exon_list = exons  # type: ignore
        else:
            for i in exon_index:
                if isinstance(i, int):
                    # exon -> blocks position
                    if i < 1 or i > len(exons):
                        raise IndexError(f"{i} is out of exon index")
                    i = exons[i - 1]
                exon_list.append(i)

        # check
        for i in exon_list:
            block = self.blocks[self._get_block_index(i)]
            if block.type != "exon":
                raise IndexError(f"{block} is not exon: input={i}")
        return self.select_block(exon_list)

    def split(self: GenemsaType) -> list[GenemsaType]:
        """Split the msa by blocks"""
        return [self.select_block([i]) for i in range(len(self.blocks))]

    # block + column function
    def shrink(self: GenemsaType) -> GenemsaType:
        """Remove empty base if all bases in that column is gap"""
        # index to delete
        freqs = self.calculate_frequency()
        masks = [f[4] != sum(f) for f in freqs]
        new_msa = self.copy(copy_allele=False)

        # recalcuate blocks
        for i in range(len(self.blocks)):
            start, end = self.get_block_interval(i)
            new_msa.blocks[i].length = sum(masks[start:end])
        assert sum(masks) == new_msa.get_length()
        new_msa.index = [new_msa.index[i] for i in range(len(masks)) if masks[i]]

        # remove base in allele
        for allele, seq in self.alleles.items():
            new_msa.alleles[allele] = "".join(
                seq[i] for i in range(len(seq)) if masks[i]
            )

        return new_msa

    def reset_index(self: GenemsaType) -> GenemsaType:
        """
        Reset index:
        The old position information will be discard.

        Each position information will be counted from 0 and
        the label and name will copy from its block information

        Example:
        ``` python
        >>> print(a.format_alignment_diff())
                           101 103 105 107 109 111 113 115 117 119
                             |   |   |   |   |   |   |   |   |   |
         A*01:01:01:01       G   G   T   C   C   A   C   C   G   A
         A*01:01:01:02N      -   -   -   -   -   -   -   -   -   -

        >>> print(a.reset_index().format_alignment_diff())
                            1
                            |
         A*01:01:01:01      GGTCCACCGA
         A*01:01:01:02N     ----------
        ```
        """
        new_msa = self.copy()
        start = 0
        new_msa.index = []
        for block in self.blocks:
            for _ in range(block.length):
                new_msa.index.append(
                    IndexInfo(
                        pos=start,
                        type=block.type,
                        name=block.name,
                    )
                )
                start += 1
        assert start == self.get_length()
        return new_msa

    def __add__(self: GenemsaType, msa: GenemsaType) -> GenemsaType:
        """
        Concat 2 MSA

        Example:
          >>> print(a_gen.select_exon([2]) + a_gen.select_exon([3]))
          <A nuc alleles=4100 block=exon2(335) exon3(413)>
        """
        names0 = set(self.get_sequence_names())
        names1 = set(msa.get_sequence_names())
        if names0 != names1:
            raise ValueError(
                "Can not concat because some allele is miss: "
                + str(names0.symmetric_difference(names1))
            )
        new_msa = self.copy()
        new_msa.blocks.extend(copy.deepcopy(msa.blocks))
        new_msa.index.extend(copy.deepcopy(msa.index))
        for name, seq in msa.alleles.items():
            new_msa.alleles[name] += seq
        return new_msa

    def __getitem__(self: GenemsaType, index: Any = None) -> GenemsaType:
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
            return self.copy()

        # Extract specific region in alignment
        Genemsa = type(self)
        if isinstance(index, int):
            index = [index]
        if isinstance(index, slice):
            new_msa = Genemsa(self.gene_name, reference=self.reference)
            new_msa.alleles = {
                allele: seq[index] for allele, seq in self.alleles.items()
            }
            new_msa.blocks = [BlockInfo(length=len(new_msa.get_reference()[1]))]
            new_msa.index = copy.deepcopy(self.index[index])
            return new_msa
        elif isinstance(index, (tuple, list)):
            new_msa = Genemsa(self.gene_name, reference=self.reference)
            new_msa.alleles = {
                allele: "".join([seq[i] for i in index])
                for allele, seq in self.alleles.items()
            }
            new_msa.blocks = [BlockInfo(length=len(new_msa.get_reference()[1]))]
            new_msa.index = copy.deepcopy([self.index[i] for i in index])
            return new_msa
        # Fail
        else:
            raise TypeError("Bad usage")

    # type-related functions
    def select_complete(self: GenemsaType) -> GenemsaType:
        """Select non exon-only sequences (No `E` in the sequence)"""
        new_msa = self.copy(copy_allele=False)
        new_msa.alleles = {
            allele: seq for allele, seq in self.alleles.items() if "E" not in seq
        }
        return new_msa

    def select_incomplete(self: GenemsaType) -> GenemsaType:
        """Select exon-only sequences (`E` exists in the sequence)"""
        new_msa = self.copy(copy_allele=False)
        new_msa.alleles = {
            allele: seq for allele, seq in self.alleles.items() if "E" in seq
        }
        return new_msa

    def fill_incomplete(self: GenemsaType, ref_allele: str) -> GenemsaType:
        """Fill the `E` in exon-only sequences with ref_allele sequence (inplace)"""
        if ref_allele not in self.alleles:
            raise KeyError(f"{ref_allele} not found")

        ref_seq = self.alleles[ref_allele]
        for allele, seq in self.alleles.items():
            if "E" in seq:
                self.alleles[allele] = "".join(
                    [seq[i] if seq[i] != "E" else ref_seq[i] for i in range(len(seq))]
                )
        return self

    def merge_exon(self: GenemsaType, msa_nuc: GenemsaType) -> GenemsaType:
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
        nuc_name_index = {
            b.name: i for i, b in enumerate(msa_nuc.blocks) if b.type == "exon"
        }

        # check it's one-to-one mapping
        exon_set = set(b.name for b in self.blocks if b.type == "exon")
        if set(nuc_name_index.keys()) != exon_set:
            raise ValueError(
                f"Cannot match blocks: " f"gen={exon_set} nuc={nuc_name_index.keys()}"
            )

        # create new msa and make sure the order of alleles
        new_msa = self.copy(copy_allele=False)
        new_msa.blocks = []
        new_msa.index = []
        new_msa.alleles = {name: "" for name in self.get_sequence_names()}
        new_msa.alleles.update({name: "" for name in msa_nuc.get_sequence_names()})

        # allele names
        gen_names = set(self.get_sequence_names())
        nuc_names = set(msa_nuc.get_sequence_names())
        exclude_name: set[str] = set()

        # block-wise
        msas_gen = self.split()
        msas_nuc = msa_nuc.split()
        for i_gen in range(len(self.blocks)):
            # intron -> fill with E
            if self.blocks[i_gen].name not in nuc_name_index:
                for name in nuc_names - gen_names:
                    msas_gen[i_gen].append(name, "E" * self.blocks[i_gen].length)
                new_msa += msas_gen[i_gen].remove(list(exclude_name))
            # exon -> check before merge
            else:
                i_nuc = nuc_name_index[self.blocks[i_gen].name]
                # content change or length change
                if msas_nuc[i_nuc].get_length() != msas_gen[i_gen].get_length() or any(
                    msas_nuc[i_nuc].get(name) != msas_gen[i_gen].get(name)
                    for name in (nuc_names & gen_names)
                ):
                    # check before merge
                    if len(gen_names - nuc_names):
                        raise ValueError(
                            f"Some alleles doesn't exist in nuc MSA: "
                            f"{gen_names - nuc_names}"
                        )

                    diff_name = filter(
                        lambda name: msas_nuc[i_nuc].get(name).replace("-", "")
                        != msas_gen[i_gen].get(name).replace("-", ""),
                        gen_names,
                    )
                    diff_names = list(diff_name)
                    if diff_names:
                        self.logger.warning(
                            f"Some exon sequences in gen MSA "
                            f"is not same as in nuc MSA "
                            f"{self.blocks[i_gen].name}: {diff_names}"
                        )
                        new_msa.remove(diff_names)
                        exclude_name.update(diff_names)
                new_msa += msas_nuc[i_nuc].remove(list(exclude_name))
        return new_msa.reset_index()

    def assume_label(self: GenemsaType, seq_type: str = "gen") -> GenemsaType:
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
            self.blocks[0].name = "5UTR"
            self.blocks[-1].type = "three_prime_UTR"
            self.blocks[-1].name = "3UTR"
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
