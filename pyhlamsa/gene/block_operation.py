from typing import List, TypeVar, Optional, Any

from .base import GenemsaBase

GenemsaType = TypeVar("GenemsaType", bound="GenemsaBlockOp")


class GenemsaBlockOp(GenemsaBase):
    """
    The class inherited from base model
    provided some block-wise operation
    """
    def get_length(self) -> int:
        """ Get the length of MSA """
        # 0 sequences is allow
        return sum(self.get_block_length())

    def get_block_length(self) -> List[int]:
        """ Get the block's length of MSA """
        return [i.length for i in self.blocks]

    def _get_block_position(self) -> List[int]:
        """ Calculate the start position of each block """
        pos = [0]
        for block in self.blocks:
            pos.append(pos[-1] + block.length)
        return pos

    def select_exon(self: GenemsaType,
                    exon_index: Optional[List[Any]] = None) -> GenemsaType:
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
        exon_list = []  # type: List[int]
        if not exon_index:
            exon_list = possible_exon_index
        else:
            for i in exon_index:
                if isinstance(i, int):
                    exon_list.append(i * 2 - 1)
                else:
                    exon_list.append([b.name for b in self.blocks].index(i))

        # check
        for ind in exon_list:
            if ind not in possible_exon_index:
                raise IndexError(f"You select the block is not exon: {ind}")
        new_msa = self.select_block(exon_list)
        return new_msa

    def select_block(self: GenemsaType, index: List[Any]) -> GenemsaType:
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
        index_ids = []
        for ind in index:
            if isinstance(ind, str):
                index_ids.append(block_name_mapping[ind])
            elif ind == -1:
                index_ids.append(len(self.blocks) - 1)
            else:
                index_ids.append(ind)

        # check index boundary
        if not (max(index_ids) < len(self.blocks) and min(index_ids) >= -1):
            raise IndexError("Check block index is correct")

        # new a msa object
        new_msa = self.copy(copy_allele=False)
        gen_pos = self._get_block_position()

        # choose block index by index
        new_block = []
        new_index = []
        for i in index_ids:
            new_block.append(self.blocks[i])
            new_index.extend(self.index[gen_pos[i]:gen_pos[i + 1]])
        new_msa.blocks = new_block
        new_msa.index = new_index

        # extract the sequences inside block region
        for allele, gen_seq in self.alleles.items():
            new_seq = "".join([gen_seq[gen_pos[i]:gen_pos[i + 1]] for i in index_ids])
            new_msa.alleles[allele] = new_seq
        return new_msa

    def split(self: GenemsaType) -> List[GenemsaType]:
        """ Split the msa by blocks """
        return [self.select_block([i]) for i in range(len(self.blocks))]
