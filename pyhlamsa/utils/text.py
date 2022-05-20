"""
Utility for generation text output string
"""
from __future__ import annotations
import dataclasses
from collections.abc import Iterator
from typing import List, TYPE_CHECKING
if TYPE_CHECKING:
    from pyhlamsa import Genemsa, IndexInfo


@dataclasses.dataclass
class _Column:
    type: str            # print the position
    length: int = 1      # type of printable word ("base" | "allele_name" | "char")
    index: bool = False  # the width
    word: str = ""       # the word for print (for type=char)
    pos: int = -1        # pos = position of base (for type=base)


def _apply_column_format_oneline(self: Genemsa, columns_format: List[_Column]) -> str:
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


def msa_to_string(self: Genemsa, **kwargs) -> str:
    """ Turn msa to string """
    output_str = ""
    for columns_format in _generate_column_format(self.index, **kwargs):
        output_str += _apply_column_format_oneline(self, columns_format)
        output_str += "\n\n"
    return output_str
