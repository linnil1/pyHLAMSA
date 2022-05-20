"""
This module is to convert msa into string

All string formatting code are written in this file
"""
import dataclasses
from collections.abc import Iterator
from typing import List, Iterable, TypeVar

from .base import IndexInfo
from .column_operation import GenemsaColumnOp


GenemsaType = TypeVar("GenemsaType", bound="GenemsaTextOp")


@dataclasses.dataclass
class _Column:
    """ Column format """
    type: str            # print the position
    length: int = 1      # type of printable word ("base" | "allele_name" | "char")
    index: bool = False  # the width
    word: str = ""       # the word for print (for type=char)
    pos: int = -1        # pos = position of base (for type=base)


def _apply_column_format_oneline(self: GenemsaType, columns_format: List[_Column]) -> str:
    """
    This function only controlling the layout not logic

    It will print the msa with the format defined in columns_format

    It may overlap the text if the length of `_Column` doesn't given enough space.
    """
    output_str = ""
    # index line
    index_left_space = 0
    for column in columns_format:
        index_left_space += column.length
        if column.index:
            # note: 1-base
            output_str += f"{self.index[column.pos].pos + 1:>{index_left_space}}"
            index_left_space = 0
    output_str += "\n"

    # indicator line
    index_left_space = 0
    for column in columns_format:
        index_left_space += column.length
        if column.index:
            output_str += f"{'|':>{index_left_space}}"
            index_left_space = 0
    output_str += "\n"

    # allele line
    for name, seq in self.alleles.items():
        for column in columns_format:
            if column.type == "char":
                output_str += column.word
            elif column.type == "allele_name":
                output_str += f"{name:<{column.length}}"
            elif column.type == "base":
                output_str += seq[column.pos]
        output_str += "\n"
    return output_str


def _generate_column_format(index: List[IndexInfo],
                            show_position_set=None,
                            wrap=100) -> Iterator[List[_Column]]:
    """ Determine the column by the index information (and block) """
    # init
    pos = 0  # current sequence position
    last_pos = -1  # last position
    last_block_name = ""  # block name of last position
    seq_length = len(index)
    column_format = []  # type: list[_Column]
    while pos < seq_length:
        show_index = False  # if this base need index

        # init a line
        if not column_format:
            column_format.append(_Column(type="char", word=" "))
            column_format.append(_Column(type="allele_name", length=18))
            column_format.append(_Column(type="char", word=" "))
            line_pos = 20  # current position in the line
            last_index_pos = 0
            num_base_in_line = 0
            show_index = True  # Should this position add position number on header

        # block indicator (but skip first block)
        if pos and last_block_name != index[pos].name:
            column_format.append(_Column(type="char", word="|"))
            line_pos += 1
            show_index = True
        last_block_name = index[pos].name

        # not consecutive position: show the index
        if index[pos].pos != last_pos + 1:
            column_format.append(_Column(type="char", word=" "))
            line_pos += 1
            show_index = True
        last_pos = index[pos].pos

        # force to show
        if show_position_set is not None:
            show_index = index[pos].pos in show_position_set

        # space to avoid index overlapping
        # pos is 0-base
        while (show_index and last_index_pos + len(str(index[pos].pos + 1)) > line_pos):
            column_format.append(_Column(type="char", word=" "))
            line_pos += 1

        # base
        column_format.append(_Column(type="base", pos=pos, index=show_index))
        pos += 1
        line_pos += 1
        num_base_in_line += 1
        if show_index:
            last_index_pos = line_pos

        # break the line
        if (num_base_in_line >= wrap or (line_pos > 120 and num_base_in_line % 10 == 0)):
            yield column_format
            column_format = []
            continue

        # split 10 bases in a line
        if num_base_in_line % 10 == 0:
            column_format.append(_Column(type="char", word=" "))
            line_pos += 1

    if column_format:
        yield column_format


def msa_to_string(self: GenemsaType, **kwargs) -> str:
    """ Turn msa to string """
    output_str = ""
    for columns_format in _generate_column_format(self.index, **kwargs):
        output_str += _apply_column_format_oneline(self, columns_format)
        output_str += "\n\n"
    return output_str


class GenemsaTextOp(GenemsaColumnOp):
    """ The class is to transfer msa into string """
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
        if isinstance(pos, int):
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
            if merged_bases and merged_bases[-1][1] + right * 2 >= b:
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
