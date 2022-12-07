import dataclasses
from typing import Any, TypeVar, Type
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment

from .base import IndexInfo, BlockInfo
from .column_operation import GenemsaColumnOp


GenemsaType = TypeVar("GenemsaType", bound="GenemsaConverter")


class GenemsaConverter(GenemsaColumnOp):
    """
    This class is to convert msa into another format and convert them back

    Support format:
    * dict(json)
    * SeqRecord
    * MultipleSeqAlignment
    """
    def meta_to_json(self) -> dict[str, Any]:
        """ Extract all meta information about this msa into json """
        meta = {
            'index': [dataclasses.asdict(i) for i in self.index],
            'blocks': [dataclasses.asdict(b) for b in self.blocks],
            'name': self.gene_name,
            'reference': self.reference,
        }
        return meta

    @classmethod
    def meta_from_json(cls: Type[GenemsaType],
                       data: dict[str, Any] = {}) -> GenemsaType:
        """ Import meta information from json """
        Genemsa = cls
        if data:
            return Genemsa(data['name'],
                           blocks=[BlockInfo(**b) for b in data['blocks']],
                           index=[IndexInfo(**i) for i in data['index']],
                           reference=data.get("reference"))
        return Genemsa("Unamed")

    def to_records(self: GenemsaType, gap=True) -> list[SeqRecord]:
        """
        Transfer MSA to list of SeqRecord

        Args:
          gap (bool): The sequence included gap or not
        """
        if gap:
            return [SeqRecord(Seq(seq.replace("E", "-")),
                              id=allele, description="")
                    for allele, seq in self.alleles.items()]
        else:
            return [SeqRecord(Seq(seq.replace("E", "").replace("-", "")),
                              id=allele, description="")
                    for allele, seq in self.alleles.items()]

    def to_MultipleSeqAlignment(self) -> MultipleSeqAlignment:
        """ Transfer this object to MultipleSeqAlignment(biopython) """
        return MultipleSeqAlignment(self.to_records(gap=True),
                                    annotations=self.meta_to_json())

    @classmethod
    def from_MultipleSeqAlignment(cls: Type[GenemsaType],
                                  bio_msa: MultipleSeqAlignment) -> GenemsaType:
        """
        Transfer MultipleSeqAlignment instance(biopython) to Genemsa

        See more details in [biopython](
        https://biopython.org/docs/1.75/api/Bio.Align.html#Bio.Align.MultipleSeqAlignment)
        """
        new_msa = cls.meta_from_json(bio_msa.annotations)
        if not new_msa.blocks:
            new_msa.blocks = [BlockInfo(length=bio_msa.get_alignment_length())]
        if not new_msa.index:
            new_msa = new_msa.reset_index()

        for seq in bio_msa:
            new_msa.alleles[seq.id] = str(seq.seq)
        assert new_msa.get_length() == bio_msa.get_alignment_length()
        return new_msa
