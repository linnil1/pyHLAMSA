"""
This module is to convert msa into string

All string formatting code are written in this file
"""
import dataclasses
from collections.abc import Iterator, Iterable
from typing import TypeVar, Union, Any

from itertools import chain

from .base import IndexInfo, GenemsaBase


GenemsaType = TypeVar("GenemsaType", bound="GenemsaTextOp")


@dataclasses.dataclass
class _TextColumn:
    """
    Column format.

    The list of this format define the layout of the page

    type:
        * char: print character (in value)
        * allele_name: print allele name with reserved length
        * index: show the index (the position in value)
        * indicator: show the index (the position in value)
    value: The main value (The meaning is defined by the type)
    length: The width that the column need
    """

    type: str  # ch
    value: Any
    length: int = 1  # length of the format


_TextPage = list[_TextColumn]


class GenemsaTextOp(GenemsaBase):
    """
    The class is to transfer msa into string.

    Those methods shares similar structures:
    1. Input MSA
    2. Modfided MSA for specific aim
    3. Calculate MSA's block/index to page format
    4. Apply MSA into page format
    5. Print it
    """

    def _format_page(
        self,
        show_position_set: set[int] = set(),
        max_page_width: int = 140,
    ) -> Iterator[_TextPage]:
        """
        Step3: Determine when to space, when to split, where to breakline.
        And return a list of _TextColumn.

        Rules:
            0. Calculate where to print sequences
            1. Start with ' {allele_name:18s} '
            2. Split every different blocks
            3. Break every 10 bases
            4. Position indicator exists when line or block start or position is not continued
            5. Position indicator can be forcely turn on or turn off, left-aligned and right-aligned
            6. Page should but over the screen
        """
        # init
        seq_index = 0  # current sequence position
        seq_length = self.get_length()  # sequence max length
        last_seq_text_value = -1  # Last position of sequence's shwoing position(text)
        # (Note the position is not seq_index, it's absolute position of origin MSA)
        last_block_name = "last_block_name"  # Block name of last base

        page: _TextPage = []
        while seq_index < seq_length:
            show_text = False  # if this base need column position indicator

            # init a line
            if not page:
                # rule 1
                page.append(_TextColumn(type="char", value=" "))
                page.append(_TextColumn(type="allele_name", value=18, length=18))
                page.append(_TextColumn(type="char", value=" "))
                column_pos = 20  # current position in the line
                last_avail_text_pos = 0  # the position allow to print column position
                num_base_in_line = (
                    0  # number of base in the line, mostly I'll split by 10
                )
                # rule4
                show_text = True  # Should this position add position number on header

            # block indicator (rule2, but skip first block)
            if seq_index and last_block_name != self.index[seq_index].name:
                page.append(_TextColumn(type="char", value="|"))
                column_pos += 1
                show_text = True
            last_block_name = self.index[seq_index].name

            # rule3
            if num_base_in_line % 10 == 0:
                page.append(_TextColumn(type="char", value=" "))
                column_pos += 1

            # rule6
            if column_pos >= max_page_width or (
                column_pos > max_page_width - 20 and num_base_in_line % 10 == 0
            ):
                yield page
                page = []
                continue

            # rule4
            if self.index[seq_index].pos != last_seq_text_value + 1:
                page.append(_TextColumn(type="char", value=" "))
                column_pos += 1
                show_text = True
            last_seq_text_value = self.index[seq_index].pos

            # rule5
            if show_position_set:
                show_text = self.index[seq_index].pos in show_position_set

            # rule5
            # space to avoid index overlapping
            # pos is 0-base and printed postion text is 1-base
            # TODO: left-align and right-align position text
            text = str(self.index[seq_index].pos + 1)
            while show_text and last_avail_text_pos + len(text) > column_pos:
                page.append(_TextColumn(type="char", value=" "))
                column_pos += 1

            # rule0, rule4
            page.append(_TextColumn(type="base", value=seq_index))
            if show_text:
                page.append(_TextColumn(type="indicator", value=True, length=0))
                # index has it's own format
                page.append(
                    _TextColumn(
                        type="index",
                        value=text,
                        length=column_pos + 1 - last_avail_text_pos,
                    )
                )
                last_avail_text_pos = column_pos + 1
            seq_index += 1
            column_pos += 1
            num_base_in_line += 1

        if page:
            # Not need to add "\n" every line in here
            # add the "\n" in _apply_page
            # page.append(_TextColumn(type="char", value="\n"))
            yield page

    def _apply_page(self, page: _TextPage) -> Iterator[str]:
        """Step 4: Page + MSA -> str"""
        output_str = ""

        # column position line
        for column in page:
            if column.type == "index":
                # note: 1-base
                output_str += f"{column.value:>{column.length}}"
        output_str += "\n"
        yield output_str

        # indicator line
        output_str = ""
        index_left_space = 0
        for column in page:
            if column.type != "index":
                index_left_space += column.length
            if column.type == "indicator":
                output_str += f"{'|':>{index_left_space}}"
                index_left_space = 0
        output_str += "\n"
        yield output_str

        # allele line
        for name, seq in self.alleles.items():
            output_str = ""
            for column in page:
                if column.type == "char":
                    output_str += column.value
                elif column.type == "allele_name":
                    output_str += f"{name:<{column.length}}"
                elif column.type == "base":
                    output_str += seq[column.value]
            output_str += "\n"
            yield output_str

        # separte each page by two lines
        yield "\n\n"

    def format_alignment(self, **kwargs_format) -> Iterator[str]:
        """
        Transfer MSA to string(generator)
        To the things in step3 and step4.

        Args:
            kwargs_format:
              * show_position_set(set[int]): Force to show the position
              * max_page_width(int): The max width of each line (Default = 140char)
        """
        if not self:
            raise ValueError("MSA is empty")
        pages = self._format_page(**kwargs_format)
        yield from chain.from_iterable(map(self._apply_page, pages))

    def print_alignment(self, **kwargs_format) -> None:
        """
        Print the MSA

        Note: The index shown on output string is 1-base absolute position.
        If you want to print in relative position, use `.reset_index()` first.

        Args:
          wrap: The maximum base per line
          show_position_set (set): The list of position(0-base) you want to indicate

        Returns:
          str: A formatted string

        Examples:
          >>> a = msa.select_allele(r"A\\*.*:01:01:01$").select_exon([6,7])
          >>> print("\n".join(a.format_alignment()))
          # or simply use
          >>> a.print_alignment()
                              3166                                  3342
                                 |                                     |
             A*01:01:01:01       ATAGAAAAGG AGGGAGTTAC ACTCAGGCTG CAA| GCAGTGA CAGTGCCCAG
             A*02:01:01:01       ATAGAAAAGG AGGGAGCTAC TCTCAGGCTG CAA| GCAGTGA CAGTGCCCAG
             A*03:01:01:01       ATAGAAAAGG AGGGAGTTAC ACTCAGGCTG CAA| GCAGTGA CAGTGCCCAG
             A*11:01:01:01       ATAGAAAAGG AGGGAGTTAC ACTCAGGCTG CAA| GCAGTGA CAGTGCCCAG
             A*23:01:01:01       ATAGAAAAGG AGGGAGCTAC TCTCAGGCTG CAA| GCAGTGA CAGTGCCCAG
             A*25:01:01:01       ATAGAAAAGG AGGGAGCTAC TCTCAGGCTG CAA| GCAGTGA CAGTGCCCAG
        """
        page_str = self.format_alignment(**kwargs_format)
        for i in page_str:
            print(i, end="")

    def print_alignment_diff(self, ref_allele: str = "", **kwargs_format) -> None:
        """
        Print the sequences of all alleles diff from `ref_allele` sequence.

        * `-` indicate same as reference
        * `*` indicate deletion
        * `ATCG` indicate the SNV
        * `E` indicate the null base in exon-only sequence
        * `|` is intron and exon boundary

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
        new_msa = self._calc_diff_msa(ref_allele)
        new_msa.print_alignment(**kwargs_format)

    def _calc_diff_msa(self: GenemsaType, ref_allele: str = "") -> GenemsaType:
        """Step 2. Modifiy the MSA. See print_alignment_diff for result"""
        ref_allele, ref_seq = self.get_allele_or_error(ref_allele)

        # use new msa object to save sequences
        new_msa = self.copy(copy_allele=False)
        new_msa.alleles = {ref_allele: ref_seq.replace("-", "*")}
        for allele, seq in self.items():
            if allele == ref_allele:
                continue
            new_seq = ""
            for i in range(len(seq)):
                if seq[i] == "-":
                    new_seq += "*"
                elif seq[i] == ref_seq[i]:
                    new_seq += "-"
                else:
                    new_seq += seq[i]
            new_msa.alleles[allele] = new_seq
        return new_msa

    def _calc_center_msa(
        self: GenemsaType, pos: Iterable[int], left: int = 5, right: int = 5
    ) -> GenemsaType:
        """Step 2. Modifiy the MSA. See print_alignment_from_center for result"""
        msa = None
        for p in want_pos:
            if msa is None:
                msa = self[p - left : p + right]
            else:
                msa += self[p - left : p + right]
        assert msa
        return msa

    def print_alignment_from_center(
        self: GenemsaType,
        pos: Union[int, Iterable[int]],
        left: int = 5,
        right: int = 5,
        **kwargs_format,
    ) -> GenemsaType:
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
        new_msa = self._calc_center_msa(want_pos, left, right)
        show_position_set = set(self.index[i].pos for i in want_pos)
        new_msa.print_alignment(show_position_set=show_position_set, **kwargs_format)

    def print_snv(self, **kwargs_format) -> None:
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
        msa, grouped_bases = self._calc_snp_msa()
        show_position_set = set(
            self.index[i].pos for bases in grouped_bases for i in bases
        )
        print(f"#Total variantion: {sum(map(len, grouped_bases))}\n")
        msa = msa._calc_diff_msa()
        msa.print_alignment(show_position_set=show_position_set, **kwargs_format)

    def _calc_snp_msa(self: GenemsaType) -> tuple[GenemsaType, list[list[int]]]:
        bases = self.get_variantion_base()
        if not bases:
            return type(GenemsaType)("Empty"), []

        # merge if two variant are too close
        grouped_bases: list[list[int]] = []
        right = 5
        for b in bases:
            if grouped_bases and grouped_bases[-1][-1] + right * 2 >= b:
                grouped_bases[-1].append(b)
            else:
                grouped_bases.append([b, b])

        # format the string
        # Using the property of index
        msa = None
        for grouped_base in grouped_bases:
            b_left, b_right = grouped_base[0], grouped_base[-1]
            if msa is None:
                msa = self[max(b_left - 5, 0) : b_right + 5]
                print("Here", msa)
                print("Here", msa.get_length())
            else:
                msa += self[max(b_left - 5, 0) : b_right + 5]
        assert msa
        return msa, grouped_bases
