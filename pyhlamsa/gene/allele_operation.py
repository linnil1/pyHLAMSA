import re
import copy
from typing import List, TypeVar, Iterable, Tuple, Union
from Bio.Seq import Seq

from ..utils import cigar
from .block_operation import GenemsaBlockOp

GenemsaType = TypeVar("GenemsaType", bound="GenemsaAlleleOp")


class GenemsaAlleleOp(GenemsaBlockOp):
    """
    This class provided method for allele-wise (row-wise) operation
    """
    def get(self, allele: str) -> str:
        """ Get the sequence by allele name """
        return self.alleles[allele]

    def sort(self: GenemsaType) -> GenemsaType:
        """ Sort the sequences """
        self.alleles = dict(sorted(self.alleles.items(), key=lambda i: i[1]))
        return self

    def sort_name(self: GenemsaType) -> GenemsaType:
        """ Sort the sequences by name """
        self.alleles = dict(sorted(self.alleles.items(), key=lambda i: i[0]))
        return self

    def append(self: GenemsaType, name: str, seq: str) -> GenemsaType:
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

    def extend(self: GenemsaType, msa: GenemsaType) -> GenemsaType:
        """ Add MSA's alleles into this MSA (inplace) """
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

    def select_allele(self: GenemsaType, query: Union[str, List[str]]) -> GenemsaType:
        """
        Select allele name by regex or list of name

        Examples:
          >>> # select allele name start with A*01:01
          >>> msa.select_allele("A\\*01\\:01.*")
          >>> # select allele by list of string
          >>> msa.select_allele(["A*01:01", "A*02:03"])
        """
        new_msa = self.copy(copy_allele=False)
        if isinstance(query, str):
            new_msa.alleles = {allele: seq
                               for allele, seq in self.alleles.items()
                               if re.match(query, allele)}
        elif isinstance(query, list):
            new_msa.alleles = {name: self.alleles[name] for name in query}
        return new_msa

    def reverse_complement(self: GenemsaType) -> GenemsaType:
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
          `cigar = [(M, 1), (X, 1), (D, 2), (M, 2)]`
        """
        if target_allele not in self.alleles:
            raise KeyError(f"{target_allele} not found")
        if not ref_allele:
            ref_allele = self.get_reference()[0]
        if ref_allele not in self.alleles:
            raise KeyError(f"{ref_allele} not found")

        return cigar.calculate_cigar(self.alleles[ref_allele],
                                     self.alleles[target_allele])
