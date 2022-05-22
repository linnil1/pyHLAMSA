import copy
from typing import List, TypeVar

from .base import BlockInfo, IndexInfo
from .block_operation import GenemsaBlockOp


GenemsaType = TypeVar("GenemsaType", bound="GenemsaColumnOp")


class GenemsaColumnOp(GenemsaBlockOp):
    """
    The class inherited from base model
    provided some column-wise and index-wise operation
    """

    def calculate_frequency(self) -> List[List[int]]:
        """
        Calculate ATCG and gap frequency of each bp in MSA

        Returns:
          frequency (list[list[int]]):
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
                    "MSA contains gap, try .shrink() before .get_consensus()")
            max_ind = [max(range(4), key=lambda i: f[i]) for f in freqs]
        else:
            max_ind = [max(range(5), key=lambda i: f[i]) for f in freqs]
        seq = ["ATCG-"[i] for i in max_ind]
        return "".join(seq)

    def shrink(self: GenemsaType) -> GenemsaType:
        """ Remove empty base if all bases in that column is gap """
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
            raise ValueError("Can not concat because some allele is miss: "
                             + str(names0.symmetric_difference(names1)))
        new_msa = self.copy()
        new_msa.blocks.extend(copy.deepcopy(msa.blocks))
        new_msa.index.extend(copy.deepcopy(msa.index))
        for name, seq in msa.alleles.items():
            new_msa.alleles[name] += seq
        return new_msa

    def __getitem__(self: GenemsaType, index=None) -> GenemsaType:
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
            new_msa.alleles = {allele: seq[index]
                               for allele, seq in self.alleles.items()}
            new_msa.blocks = [BlockInfo(length=len(new_msa.get_reference()[1]))]
            new_msa.index = copy.deepcopy(self.index[index])
            return new_msa
        elif isinstance(index, (tuple, list)):
            new_msa = Genemsa(self.gene_name, reference=self.reference)
            new_msa.alleles = {allele: "".join([seq[i] for i in index])
                               for allele, seq in self.alleles.items()}
            new_msa.blocks = [BlockInfo(length=len(new_msa.get_reference()[1]))]
            new_msa.index = copy.deepcopy([self.index[i] for i in index])
            return new_msa
        # Fail
        else:
            raise TypeError("Bad usage")

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
                new_msa.index.append(IndexInfo(
                    pos=start,
                    type=block.type,
                    name=block.name,
                ))
                start += 1
        assert start == self.get_length()
        return new_msa
