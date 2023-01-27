from __future__ import annotations
from typing import TypeVar, Union, Any, TextIO, Type
from itertools import chain
import os
import json
import logging
import dataclasses

from Bio import AlignIO, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
import pysam
from pysam import bcftools

from .base import GenemsaBase, BlockInfo, IndexInfo
from .text_operation import GenemsaTextOp
from ..utils import vcf, alignment, cigar


logger = logging.getLogger(__name__)
GenemsaType = TypeVar("GenemsaType", bound="GenemsaIO")
FileType = Union[str, os.PathLike, TextIO]


def getFileHandle(file: FileType) -> tuple[TextIO, bool]:
    """
    Return file handle from filename
    and the file object should be closed or not.
    """
    if isinstance(file, (str, os.PathLike)):
        return open(file, "w"), True
    else:
        return file, False


def _table_to_string(data: list[list[Any]]) -> str:
    """Turn table into tab separted string"""
    return "\n".join("\t".join(map(str, items)) for items in data)


def _block_to_gff(
    blocks: list[BlockInfo],
    allele: str,
    strand: str = "+",
    igv_show_label: bool = False,
) -> list[Any]:
    # Gene
    records = [
        [
            allele,
            "pyHLAMSA",
            "gene",
            1,
            sum(i.length for i in blocks),
            ".",
            strand,
            ".",
            f"ID={allele};Name={allele}",
        ]
    ]

    # Blocks
    pos = 0
    for b in blocks:
        # http://gmod.org/wiki/GFF3
        # gff3 format:
        #   1. header: ref source type start end . strand . tags
        #   2. pos: 1-base included position
        #   3. type: In HLA annotations exon=CDS
        records.append(
            [
                allele,
                "pyHLAMSA",
                b.type if b.type != "exon" else "CDS",
                pos + 1,
                pos + b.length,
                ".",
                strand,
                ".",
                f"ID={b.name}_{allele}",
            ]
        )
        # To show the label of all block in IGV
        # I break the relation (Remove parent attribute)
        if not igv_show_label:
            records[-1][-1] += f";Parent={allele}"  # type: ignore
        pos += b.length
    return records


class GenemsaIO(GenemsaTextOp, GenemsaBase):
    """
    The IO operation of Genemsa.

    The file type inclues
    * MSA  -> list of SeqRecord
    * MSA  -> fasta (fa)
    * MSA <-> dict(json) (index-only)
    * MSA <-> MultipleSeqAlignment
    * MSF  -> MSA
    * MSA  -> bam
    * MSA  -> gff
    * MSA <-> our pyHLAMSA object(fa + json)
    """

    def to_records(self, gap: bool = True) -> list[SeqRecord]:
        """
        Transfer MSA to list of SeqRecord.

        If the base is `E`
        which indicate the bases are not existed for exon-only sequence,
        those beses are treated as gap.

        Args:
          gap (bool): The sequence included gap or not
        """
        if gap:
            return [
                SeqRecord(Seq(seq.replace("E", "-")), id=allele, description="")
                for allele, seq in self.items()
            ]
        else:
            return [
                SeqRecord(
                    Seq(seq.replace("E", "").replace("-", "")),
                    id=allele,
                    description="",
                )
                for allele, seq in self.items()
            ]

    def to_fasta(
        self, fname: FileType, gap: bool = True, ref_only: bool = False
    ) -> None:
        """
        Save the MSA into fasta

        Args:
            gap (bool): The sequence included gap or not
            ref_only (bool): Save reference sequence only
        """
        if ref_only:
            msa = self.select_allele([self.get_reference()[0]])
        else:
            msa = self
        SeqIO.write(msa.to_records(gap=gap), fname, "fasta")

    def meta_to_json(self) -> dict[str, Any]:
        """Extract all meta information about this msa into dict(json)"""
        meta = {
            "index": [dataclasses.asdict(i) for i in self.index],
            "blocks": [dataclasses.asdict(b) for b in self.blocks],
            "name": self.gene_name,
            "reference": self.reference,
        }
        return meta

    @classmethod
    def meta_from_json(
        cls: Type[GenemsaType], data: dict[str, Any] = {}
    ) -> GenemsaType:
        """Import meta information from json"""
        Genemsa = cls
        if data:
            return Genemsa(
                data["name"],
                blocks=[BlockInfo(**b) for b in data["blocks"]],
                index=[IndexInfo(**i) for i in data["index"]],
                reference=data.get("reference", ""),
            )
        return Genemsa("Unamed")

    def to_MultipleSeqAlignment(self) -> MultipleSeqAlignment:
        """Transfer this object to MultipleSeqAlignment(biopython)"""
        return MultipleSeqAlignment(
            self.to_records(gap=True), annotations=self.meta_to_json()
        )

    @classmethod
    def from_MultipleSeqAlignment(
        cls: Type[GenemsaType], bio_msa: MultipleSeqAlignment
    ) -> GenemsaType:
        """
        Transfer MultipleSeqAlignment instance(biopython) to Genemsa

        See more details in [biopython](
        https://biopython.org/docs/1.75/api/Bio.Align.html#Bio.Align.MultipleSeqAlignment)
        """
        new_msa = cls.meta_from_json(bio_msa.annotations)
        # Fix header (block and index)
        if not new_msa.blocks:
            new_msa.blocks = [BlockInfo(length=bio_msa.get_alignment_length())]
        if not new_msa.index:
            new_msa = new_msa.reset_index()
        # Add sequences
        for seq in bio_msa:
            new_msa.append(seq.id, str(seq.seq))
        assert new_msa.get_length() == bio_msa.get_alignment_length()
        return new_msa

    @classmethod
    def read_msf_file(cls: Type[GenemsaType], file_msf: str) -> GenemsaType:
        """Read .msf file"""
        return cls.from_MultipleSeqAlignment(AlignIO.read(file_msf, "msf"))

    def to_bam(self, fname: str, ref_allele: str = "", save_ref: bool = True) -> None:
        """
        Save the MSA into bam

        All alleles will seen as reads aligned on `ref_allele`

        Args:
          fname (str): The name of bamfile
          ref_allele (str): The reference allele.
              If the ref_allele is empty, the first allele will be reference.
          save_ref (bool): The reference allele will also be saved in the bam file
        """
        ref_allele, ref_seq = self.get_allele_or_error(ref_allele)

        # setup reference and header
        header = {
            "HD": {"VN": "1.0"},
            "SQ": [
                {"LN": len(ref_seq.replace("-", "").replace("E", "")), "SN": ref_allele}
            ],
        }

        # write bam file
        with pysam.AlignmentFile(fname, "wb", header=header) as outf:
            for allele, seq in self.items():
                # skip
                if not save_ref and allele == ref_allele:
                    continue

                # init bam record
                a = pysam.AlignedSegment()
                a.query_name = allele
                a.query_sequence = seq.replace("E", "").replace("-", "")
                a.cigar = cigar.cigar_to_pysam(self.get_cigar(allele, ref_allele))  # type: ignore

                # set to default
                a.mapping_quality = 60
                a.reference_id = 0
                a.reference_start = 0
                # a.template_length = 0
                # a.query_qualities = [30] * len(a.query_sequence)
                # a.flag = 0
                # a.next_reference_id = 0
                # a.next_reference_start = 0
                outf.write(a)

        pysam.sort("-o", fname, fname)  # type: ignore
        pysam.index(fname)  # type: ignore

    def to_gff(
        self,
        fname: FileType,
        strand: str = "+",
        ref_allele: str = "",
        igv_show_label: bool = False,
        save_all: bool = False,
    ) -> None:
        """
        Save to GFF3 format

        Args:
          fname (str): The file name of gff3
          strand (str): Must be "+" or "-".

              If the strand is "-", it will add `strand` in GFF file,
              if you want to reverse the sequence, please use `reverse_complement` first.

          ref_allele (str): The name of allele (Must be the same in save_bam)
          igv_show_label (bool): If it's false, it will generate proper GFF3.
              Set it for True as default for easiler label reading in IGV.
          save_all (bool):
              Set it True if you want to create gff records for all alleles in msa
              Note that this is not very fast.
        """
        # TODO: should I save strand == '-' in model?
        ref_allele = self.get_allele_or_error(ref_allele)[0]

        # labels
        if not all(b.type for b in self.blocks):
            self.logger.warning(
                "You should assign block's label. (We assume seq_type='other')"
            )
            self.assume_label("other")

        if save_all:
            alleles = list(self.list_alleles())
        else:
            alleles = [ref_allele]

        records = []
        for allele in alleles:
            self.logger.debug(f"{allele} to gff")

            # remove gap
            msa = self.select_allele([allele]).shrink()
            records.extend(
                _block_to_gff(
                    msa.blocks, allele, strand=strand, igv_show_label=igv_show_label
                )
            )

        # save
        f_gff, require_close = getFileHandle(fname)
        f_gff.write("##gff-version 3\n")
        f_gff.write(_table_to_string(records))
        if require_close:
            f_gff.close()

    def save_msa(self, file_fasta: FileType, file_json: FileType) -> None:
        """Save Genemsa to pyhlamsa format (fasta and json)"""
        self.to_fasta(file_fasta, gap=True)
        f_json, require_close = getFileHandle(file_json)
        json.dump(self.meta_to_json(), f_json)
        if require_close:
            f_json.close()

    @classmethod
    def load_msa(
        cls: Type[GenemsaType], file_fasta: str, file_json: str
    ) -> GenemsaType:
        """
        load Genemsa from fasta and json

        Example:
          >>> a_gen.save_msa("a_gen.fa", "a_gen.json")
          >>> a_gen = Genemsa.load_msa("a_gen.fa", "a_gen.json")
        """
        # read
        msa = AlignIO.read(file_fasta, "fasta")
        with open(file_json) as f:
            msa.annotations = json.load(f)
        # main
        new_msa = cls.from_MultipleSeqAlignment(msa)
        assert len(new_msa.get_reference()[1]) == new_msa.get_length()
        return new_msa

    def to_vcf(
        self,
        file_vcf: str,
        ref_allele: str = "",
        save_ref: bool = True,
        plain_text: bool = False,
    ) -> None:
        """
        Save Genemsa into vcf format

        It will save msa into sorted and normalized vcf file (xxx.vcf.gz)
        along with vcf's index (xxx.vcf.gz.tbi)
        (You can still output plain vcf without any manipulation by set plain_text=True)

        Note that vcf-format discard the per-base alignment especially
        when there are many snps/indels mixed together in that position.
        Also, after vcf is normalized, the msa may not be the same MSA alignment as before.
        (But sequences are still the same)

        Args:
          ref_allele (str):
            The name of reference allele.
            I recommend the reference allele contains no gaps.
          save_ref (bool): The reference allele will also be saved in the vcf file
          plain_text (bool): Disable sort vcf and index vcf
        """
        ref_allele, ref_seq = self.get_allele_or_error(ref_allele)

        # extract all variants from all alleles
        allele_variants = {}
        for allele in self.list_alleles():
            if not save_ref and allele == ref_allele:
                continue
            variants = vcf.extract_variants(ref_seq, self.get(allele))
            for v in variants:
                v.chrom = ref_allele
            allele_variants[allele] = variants

        # remove suffix
        if file_vcf.endswith(".vcf.gz"):
            basename = file_vcf[:-7]
        elif file_vcf.endswith(".vcf"):
            basename = file_vcf[:-4]
        else:
            basename = file_vcf
        tmpname = f"{basename}.tmp"

        # to vcf
        with open(f"{tmpname}.vcf", "w") as f_vcf:
            f_vcf.write(vcf.get_vcf_header(ref_allele, ref_seq))
            f_vcf.write("\n")
            f_vcf.write(_table_to_string(vcf.variants_to_table(allele_variants)))
        if plain_text:
            os.rename(f"{tmpname}.vcf", f"{basename}.vcf")
            return

        # sort, normalize and index
        self.select_allele([ref_allele]).to_fasta(f"{tmpname}.fa", gap=False)
        with open(f"{tmpname}.norm.vcf.gz", "wb") as f_vcf:
            f_vcf.write(
                bcftools.norm(  # type: ignore
                    "-f", f"{tmpname}.fa", f"{tmpname}.vcf", "-O", "z"
                )
            )
        with open(f"{basename}.vcf.gz", "wb") as f_vcf:
            f_vcf.write(bcftools.sort(f"{tmpname}.norm.vcf.gz", "-O", "z"))  # type: ignore
        bcftools.index(f"{basename}.vcf.gz", "-t", "-f")  # type: ignore
        os.remove(f"{tmpname}.fa")
        os.remove(f"{tmpname}.fa.fai")
        os.remove(f"{tmpname}.vcf")
        os.remove(f"{tmpname}.norm.vcf.gz")

    @classmethod
    def read_alignment_txt(
        cls: Type[GenemsaType], fname: str, seq_type: str = ""
    ) -> GenemsaType:
        """Same as read_imgt_alignment"""
        return cls.read_imgt_alignment(fname, seq_type)

    @classmethod
    def read_imgt_alignment(
        cls: Type[GenemsaType], fname: str, seq_type: str = ""
    ) -> GenemsaType:
        """
        Read IMGT-alignment MSA file into Genemsa object.
        e.g. IMGT/alignments/A_gen.txt
        """
        # Read all alleles
        alleles = alignment.parse_alignment_txt(fname)
        ref_seq = next(iter(alleles.values()))

        Genemsa = cls
        new_msa = Genemsa("Unnamed")
        # Use first sequence as reference
        new_msa.set_blocks([len(seq) for seq in ref_seq.split("|")])
        new_msa.assume_label(seq_type)
        new_msa.extend({name: seq.replace("|", "") for name, seq in alleles.items()})
        return new_msa

    def to_imgt_alignment(
        self,
        file: FileType,
        seq_type: str = "gen",
        index_start_from: str = "exon1",
    ) -> None:
        """
        Export to IMGT-alignment-like format.

        But some features do not implemented:
        * AA codon position: So the values in the header always 1
        * position: So the values in the header always 1

        Args:
            seq_type: gen or nuc format
            index_start_from: set the start position (1) starts from.
        """
        # This function should be located in io_operation
        new_msa = self._calc_imgt_alignment_msa()
        start = new_msa.get_block_interval(index_start_from)[0]

        # The position is calculated by reference sequence not MSA itself
        ref_no_gap_pos = -1
        seq = new_msa.get_reference()[1]
        for i, ind in enumerate(new_msa.index):
            if seq[i] != ".":
                ref_no_gap_pos += 1
            ind.pos = ref_no_gap_pos

        if seq_type == "gen":
            # force modify the index
            # IMGT ends with -1 and start with 1
            for ind in new_msa.index:
                ind.pos -= start - 1
                if ind.pos < 0:
                    ind.pos -= 1
            pages = self._format_page(
                [new_msa.index],
                header_text_align="left",
                index_header=["gDNA"],
                show_position_when_same_value=False,
            )
        elif seq_type == "nuc":
            fake_index = [IndexInfo(pos=0)] * len(new_msa.index)
            pages = self._format_page(
                [new_msa.index, fake_index],
                header_text_align="left",
                index_header=["cDNA", "AA Codon"],
                show_position_when_same_value=False,
            )
        else:
            raise NotImplementedError(f"Not implement {seq_type=}")

        # write
        txt, require_close = getFileHandle(file)
        txt.writelines(chain.from_iterable(map(new_msa._apply_page, pages)))
        if require_close:
            txt.close()
