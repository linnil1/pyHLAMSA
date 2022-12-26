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
        * index:
            show the header text (saved in value)
            If there exists more than one header, it will be index0, index1...
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

    @classmethod
    def _format_page(
        cls,
        index: list[list[IndexInfo]],
        index_header: list[str] = [],
        show_position_set: set[int] = set(),
        show_position_when_same_value: bool = True,
        position_base: int = 1,
        max_page_width: int = 140,
        header_text_align: str = "right",
        allele_name_size: int = 17,
        break_per_k_base: int = 10,
    ) -> Iterator[_TextPage]:
        """
        Step3: Determine when to space, when to split, where to breakline
        and when to showing the position text above the base.

        The function return a list of _TextColumn as a page layout.

        The layout should be like this:
        ```
         DNA(index_header)      100(index_position)  150
         RNA(index_header)        1
                                  |(indicator)         |
          A*02 (allele_name)  AATTATACCTACGGGGGAAATTTCCC(base)
          A*03 (allele_name)  AATTATACCTACGGGGGAAATTTCCC
        ```

        Rules:
            0. Calculate where to print sequences
            1. Start with ' {allele_name} '
            2. Split every different blocks
            3. Break every k bases
            4. Position text and it's indicator exists when line or block start or position is not continued
            5. index showing can be forcely turn on or turn off, left-aligned and right-aligned
            6. Page should but over the screen
            7. Index name can be written in header (e.g. gDNA)
            8. One or more index can be used and the first index is the main index
        """
        # init
        if len(index) == 0:
            raise ValueError("At least one index should be input")
        index_num = len(index)
        seq_index = 0  # current sequence position
        seq_length = len(index[0])  # sequence max length
        last_position_value = -1  # the position of the last base
        last_block_name = "last_block_name"  # block name of the last base
        column_pos = -1  # current position in the line
        last_avail_header_pos = [
            0
        ] * index_num  # the position allow to text of each index
        num_base_in_line = 0  # current number of base in the line
        last_indicator_pos = 0  # the position of last indicator

        page: _TextPage = []
        while seq_index < seq_length:
            show_index = False  # if header text should popup, so does indicator

            # init a line
            if not page:
                # rule 1
                page.append(_TextColumn(type="char", value=" "))
                page.append(
                    _TextColumn(
                        type="allele_name",
                        value=allele_name_size,
                        length=allele_name_size,
                    )
                )
                page.append(_TextColumn(type="char", value=" "))
                column_pos = 2 + allele_name_size
                last_avail_header_pos = [0] * index_num
                num_base_in_line = 0
                last_indicator_pos = 0
                show_index = True  # rule4
                if index_header:  # rule 7, 8
                    if len(index_header) != index_num:
                        raise ValueError(
                            "Number of Header is mismatched to Number of index"
                        )
                    for i in range(index_num):
                        page.append(
                            _TextColumn(
                                type=f"index{i}",
                                value=" " + index_header[i],
                                length=1 + len(index_header[i]),
                            )
                        )
                        last_avail_header_pos[i] = 1 + len(index_header[i])

            # rule2: block indicator, but skip first block
            if seq_index and last_block_name != index[0][seq_index].name:
                page.append(_TextColumn(type="char", value="|"))
                column_pos += 1
                show_index = True  # If it encounter line break, it's fine
            last_block_name = index[0][seq_index].name

            # rule3: spacing but skip split in the first base
            if num_base_in_line and num_base_in_line % break_per_k_base == 0:
                page.append(_TextColumn(type="char", value=" "))
                column_pos += 1

            # rule6: page break
            if column_pos >= max_page_width or (
                column_pos > max_page_width - 20
                and num_base_in_line % break_per_k_base == 0
            ):
                yield page
                page = []
                continue

            # rule4: position incontinuous
            if index[0][seq_index].pos != last_position_value + 1 and (
                show_position_when_same_value
                or index[0][seq_index].pos != last_position_value
            ):
                show_index = True
            last_position_value = index[0][seq_index].pos
            # rule5: Force the print position
            if show_position_set:
                show_index = index[0][seq_index].pos in show_position_set

            # rule8: for each index
            if show_index:
                for i in range(index_num):
                    text = str(index[i][seq_index].pos + position_base)
                    # rule5
                    # space to avoid index overlapping
                    while (
                        header_text_align == "right"
                        and last_avail_header_pos[i] + len(text) > column_pos
                    ) or (
                        header_text_align == "left"
                        and last_avail_header_pos[i] > column_pos - 1
                    ):
                        page.append(_TextColumn(type="char", value=" "))
                        column_pos += 1

                # rule4: Write index position
                # index has it's own format (It has independent width)
                for i in range(index_num):
                    text = str(index[i][seq_index].pos + position_base)
                    page.append(
                        _TextColumn(
                            type=f"index{i}",
                            value=text,
                            length=column_pos + 1 - last_avail_header_pos[i],
                        )
                    )
                    if header_text_align == "left":
                        page[-1].length += len(text) - 1
                        last_avail_header_pos[i] = column_pos + len(text)
                    else:
                        last_avail_header_pos[i] = column_pos + 1

                # rule4: Write indicator
                page.append(
                    _TextColumn(
                        type="indicator",
                        value="|",
                        length=column_pos + 1 - last_indicator_pos,
                    )
                )
                last_indicator_pos = column_pos + 1

            # rule0: Write sequence
            page.append(_TextColumn(type="base", value=seq_index))
            seq_index += 1
            column_pos += 1
            num_base_in_line += 1

        if page:
            # Not need to add "\n" every line in here
            # add the "\n" in _apply_page
            # page.append(_TextColumn(type="char", value="\n"))
            yield page

    def _apply_page(self, page: _TextPage) -> Iterator[str]:
        """Step 4: Page format + MSA -> str"""

        # index line (header)
        # rule8: support mutlitple index
        index_ids = sorted(
            set(column.type for column in page if column.type.startswith("index"))
        ) + ["indicator"]
        for index_id in index_ids:
            output_str = ""
            for column in page:
                if column.type == index_id:
                    output_str += f"{column.value:>{column.length}}"
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

    def format_alignment(self, **kwargs_format: Any) -> Iterator[str]:
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
        pages = self._format_page([self.index], **kwargs_format)
        yield from chain.from_iterable(map(self._apply_page, pages))

    def print_alignment(self, **kwargs_format: Any) -> None:
        """
        Print the MSA

        Note: The index shown on output string is 1-base absolute position.
        If you want to print in relative position, use `.reset_index()` first.

        * `*` indicate deletion
        * `ATCG` indicate the SNV
        * `E` indicate the null base in exon-only sequence
        * `|` is intron and exon boundary

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

    def print_alignment_diff(self, ref_allele: str = "", **kwargs_format: Any) -> None:
        """
        Print the sequences of all alleles diff from `ref_allele` sequence.

        Same as print_alignment, but
        once the base of alleles is same as the base in reference sequence,
        `-` will be used

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
        """
        Step 2. Modifiy the MSA.
        See print_alignment_diff for result
        """
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
        """
        Step 2. Modifiy the MSA.
        See print_alignment_from_center for result
        """
        msa = None
        for p in pos:
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
        **kwargs_format: Any,
    ) -> None:
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
            want_pos = list(pos)
        if not want_pos:
            return
        new_msa = self._calc_center_msa(want_pos, left, right)
        show_position_set = set(self.index[i].pos for i in want_pos)
        new_msa.print_alignment(show_position_set=show_position_set, **kwargs_format)

    def print_snv(self, **kwargs_format: Any) -> None:
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
        """
        Step 2. Modifiy the MSA.
        Find the snp in MSA and crop each block

        Return:
            msa: The modified MSA
            grouped_bases: A list of Grouped SNPs (The grouped when they are closed)
        """
        bases = self.get_variantion_base()
        if not bases:
            return type(self)("Empty"), []

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
            else:
                msa += self[max(b_left - 5, 0) : b_right + 5]
        assert msa
        return msa, grouped_bases

    def _calc_imgt_alignment_msa(
        self: GenemsaType, ref_allele: str = ""
    ) -> GenemsaType:
        """
        Step 2. Modifiy the MSA.
        Transform our MSA to IMGT-alignment-like characters.

        Used by `to_imgt_alignment`.
        """
        ref_allele, ref_seq = self.get_allele_or_error(ref_allele)

        # use new msa object to save sequences
        new_msa = self.copy(copy_allele=False)
        new_msa.alleles = {ref_allele: ref_seq.replace("-", ".")}
        for allele, seq in self.items():
            if allele == ref_allele:
                continue
            seq_list = list(seq)
            new_seq = ""
            for i in range(len(seq)):
                if seq[i] == "-" and ref_seq[i] == "-":
                    seq_list[i] = "."
                elif seq[i] == "-":
                    seq_list[i] = "*"
                elif seq[i] == ref_seq[i]:
                    seq_list[i] = "-"
            new_msa.alleles[allele] = "".join(seq_list)
        return new_msa
