"""
Define the basic class of Genemsa,
including GenemsaBase, BlockInfo, IndexInfo
"""
from typing import Tuple, List, Dict, TypeVar, Optional
import copy
import logging
import dataclasses


GenemsaType = TypeVar("GenemsaType", bound="GenemsaBase")


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


class GenemsaBase:
    """
    The basic definition of MSA

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
    def __init__(self, gene_name: str,
                 blocks: Optional[List[BlockInfo]] = None,
                 index: Optional[List[IndexInfo]] = None,
                 reference=None):
        self.gene_name = gene_name
        self.alleles = {}  # type: Dict[str, str]
        self.blocks = copy.deepcopy(blocks or [])  # intron exon length
        self.index = copy.deepcopy(index or [])
        self.logger = logging.getLogger(__name__)
        self.reference = reference

    # Show the MSA attribute
    def __repr__(self) -> str:
        block_info = " ".join([f"{b.name}({b.length})" for b in self.blocks])
        return f"<{self.gene_name} "\
               f"alleles={len(self.alleles)} "\
               f"block={block_info}>"

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

    def get_length(self) -> int:
        """ Get the length of MSA """
        # 0 sequences is allow
        return len(self.get_reference()[1])

    def list_alleles(self) -> List[str]:
        """
        List all the allele's sequence name in MSA

        Example:
           >>> a_gen.list_alleles()[:3]
           ['A*01:01:01:01', 'A*01:01:01:02N', 'A*01:01:01:03']
        """
        return list(self.alleles.keys())

    def get_sequence_names(self) -> List[str]:
        """ Same as list_allele_names """
        return self.list_alleles()

    def copy(self: GenemsaType, copy_allele=True) -> GenemsaType:
        """
        Clone a new MSA.

        Args:
          copy_allele (bool): copy the sequences
        """
        Genemsa = type(self)  # Child Type
        new_msa = Genemsa(self.gene_name,
                          blocks=self.blocks, index=self.index,
                          reference=self.reference)
        if copy_allele:
            new_msa.alleles = dict(self.alleles.items())
        return new_msa

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
        if not self.alleles:
            raise ValueError("MSA is empty")
        return next(iter(self.alleles.items()))
