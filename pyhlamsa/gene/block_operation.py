from typing import List, TypeVar, Optional, Any, Union, Tuple
from collections.abc import Iterable

from .base import GenemsaBase, BlockInfo

GenemsaType = TypeVar("GenemsaType", bound="GenemsaBlockOp")
BlockInput = Union[int, str, BlockInfo]


class GenemsaBlockOp(GenemsaBase):
    """
    The class inherited from base model
    provided some block-wise operation
    """
    def get_length(self) -> int:
        """ Get the length of MSA """
        # 0 sequences is allow
        # overwrite the original len(seq) method
        return sum(i.length for i in self.blocks)

    def _get_block_index(self, block=BlockInput) -> int:
        """ Find the index of the block """
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
        elif not isinstance(block, int):
            raise NotImplementedError(f"Type of {block} not work now")

        id = block
        if id < 0:
            id = len(self.blocks) + id
        if not (0 <= id < len(self.blocks)):
            raise IndexError(f"{block} is out of index")
        return id

    def get_block_interval(self, block: BlockInput) -> Tuple[int, int]:
        """ Calculate the start(included) and end index (excluded) of the block """
        index = self._get_block_index(block)
        start = sum(self.blocks[i].length for i in range(index))
        return start, start + self.blocks[index].length

    def get_block_position(self, block: BlockInput) -> int:
        """ Calculate the start position of the block """
        index = self._get_block_index(block)
        return sum(self.blocks[i].length for i in range(index))

    def select_exon(self: GenemsaType,
                    exon_index: Optional[Iterable[BlockInput]] = None) -> GenemsaType:
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
        exon_list = []  # type: List[BlockInput]
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
        index_ids = [self._get_block_index(i) for i in index]

        # new a msa object
        new_msa = self.copy(copy_allele=False)

        # choose block index by index
        new_block = []
        new_index = []
        all_pos = []
        for i in index_ids:
            new_block.append(self.blocks[i])
            start, end = self.get_block_interval(i)
            all_pos.append((start, end))
            new_index.extend(self.index[start: end])
        new_msa.blocks = new_block
        new_msa.index = new_index

        # extract the sequences inside block region
        for allele, gen_seq in self.alleles.items():
            new_seq = "".join([gen_seq[start:end] for start, end in all_pos])
            new_msa.alleles[allele] = new_seq
        return new_msa

    def split(self: GenemsaType) -> List[GenemsaType]:
        """ Split the msa by blocks """
        return [self.select_block([i]) for i in range(len(self.blocks))]
