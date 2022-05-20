"""
msaio: IO operation of Genemsa

Main objective
* Read msa from file(txt, msf, fa) into Genemsa object or
* Write Genemsa into file (gff, bam, fa)

usage:
``` python
from pyhlamsa import msaio
```
"""
from __future__ import annotations
import re
import json
import logging
from typing import List, Tuple, Dict
from Bio import AlignIO, SeqIO
import pysam

from .gene import Genemsa, BlockInfo


logger = logging.getLogger(__name__)


def read_alignment_txt(fname: str, seq_type="") -> Genemsa:
    """
    Read MSA file defined in IMGT .e.g. IMGT/alignments/A_gen.txt
    """
    # Read all alleles
    alleles = _parse_alignment_txt(fname)
    new_msa = Genemsa("Unnamed")
    new_msa.alleles = {name: seq.replace("|", "") for name, seq in alleles.items()}

    # Use first sequence as reference
    ref_seq = list(alleles.values())[0]
    new_msa.blocks = [BlockInfo(length=len(seq)) for seq in ref_seq.split("|")]
    new_msa.assume_label(seq_type)
    return new_msa


def _parse_alignment_txt(fname: str) -> Dict[str, str]:
    """
    Read MSA file defined in IMGT .e.g. IMGT/alignments/A_gen.txt

    Returns:
        alleles (dict of str,str): The dictionary that map allele name to sequence
    """
    # parse aligments
    alleles = {}
    ref_allele = ""
    with open(fname) as f_txt:
        for line in f_txt:
            line = line.strip()
            if not re.match(r"^\w+\*", line):
                continue

            match = re.findall(r"^(.*?) +(.*)", line)
            # assert match
            # emtpy seq (B_nuc.txt)
            if not match:
                continue
            allele, seq = match[0]

            if allele not in alleles:
                if not ref_allele:
                    ref_allele = allele
                alleles[allele] = ""

            alleles[allele] += seq

    # check sequences and replace
    rm_allele = []
    for allele in alleles:
        alleles[allele] = alleles[allele].replace(" ", "").replace("*", ".")
        if allele == ref_allele:
            alleles[allele] = alleles[allele].replace(".", "-")
            continue
        ref_seq = alleles[ref_allele]
        if len(alleles[allele]) != len(ref_seq):
            rm_allele.append(allele)
            continue
        # assert len(alleles[allele]) == len(ref_seq)

        seq = list(alleles[allele])
        for i in range(len(seq)):
            if seq[i] == "-":
                seq[i] = ref_seq[i]

            if seq[i] == "|":
                assert seq[i] == ref_seq[i]
            else:
                assert seq[i] in "ATCG."
        alleles[allele] = "".join(seq).replace(".", "-")

    for allele in rm_allele:
        logger.warning(f"{allele} violate msa length in {fname}. Removed")
        del alleles[allele]
    return alleles


def read_msf_file(file_msf: str) -> Genemsa:
    """ Read .msf file """
    return Genemsa.from_MultipleSeqAlignment(AlignIO.read(file_msf, "msf"))


def _cigar_to_pysam(cigar: List[Tuple[str, int]]) -> List[Tuple[int, int]]:
    """
    Translate cigar to cigar tuple defined in
    https://pysam.readthedocs.io/en/latest/api.html#pysam.AlignedSegment.cigartuples
    """
    op_type_map = {
        "M": 0,
        "I": 1,
        "D": 2,
        "X": 8,
    }
    return list(map(lambda i: (op_type_map.get(i[0], 0), i[1]), cigar))


def to_fasta(self: Genemsa, fname: str, gap=True):
    """
    Save the MSA into fasta

    Args:
        gap (bool): The sequence included gap or not
    """
    SeqIO.write(self.to_records(gap=gap), fname, "fasta")


def to_bam(self: Genemsa, fname: str, ref_allele="", save_ref=False):
    """
    Save the MSA into bam

    All alleles will seen as reads aligned on `ref_allele`

    Args:
      fname (str): The name of bamfile
      ref_allele (str): The reference allele.
          If the ref_allele is empty, the first allele will be reference.
      save_ref (bool): The reference allele will also be saved in the bam file
    """
    if not len(self.alleles):
        raise ValueError("MSA is empty")
    if not ref_allele:
        ref_allele = self.get_reference()[0]
    if ref_allele not in self.alleles:
        raise ValueError(f"{ref_allele} not found")
    if not fname:
        raise ValueError("filename is required")

    # setup reference and header
    ref_seq = self.alleles[ref_allele]
    header = {'HD': {'VN': "1.0"},
              'SQ': [{'LN': len(ref_seq.replace("-", "").replace("E", "")),
                      'SN': ref_allele}]}

    # write bam file
    with pysam.AlignmentFile(fname, "wb", header=header) as outf:
        for allele, seq in self.alleles.items():
            # skip
            if not save_ref and allele == ref_allele:
                continue

            # init bam record
            a = pysam.AlignedSegment()
            a.query_name = allele
            a.query_sequence = seq.replace("E", "").replace("-", "")
            a.cigar = _cigar_to_pysam(self.get_cigar(allele, ref_allele))  # type: ignore

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


def to_gff(self: Genemsa, fname: str, strand="+", ref_allele="", igv_show_label=False):
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

    """
    # TODO: should I save strand == '-' in model?
    if not len(self.blocks):
        raise ValueError("MSA is empty")
    if not ref_allele:
        ref_allele = self.get_reference()[0]
    if ref_allele not in self.alleles:
        raise ValueError(f"{ref_allele} not found")

    # labels
    if not all(b.type for b in self.blocks):
        self.logger.warning(
            "You should assign block's label. (We assume seq_type='other')")
        self.assume_label("other")

    # Gene
    gene_name = self.gene_name or ref_allele
    records = [[
        self.gene_name or ref_allele, "pyHLAMSA", "gene",
        str(1), str(self.get_length()), ".", strand, ".",
        f"ID={gene_name};Name={gene_name}"
    ]]

    # Blocks
    pos = 0
    for b in self.blocks:
        # http://gmod.org/wiki/GFF3
        # gff3 format:
        #   1. header: ref source type start end . strand . tags
        #   2. pos: 1-base included position
        #   3. type: In HLA annotations exon=CDS
        records.append(
            [ref_allele, "pyHLAMSA",
             b.type if b.type != "exon" else "CDS",
             str(pos + 1), str(pos + b.length), ".", strand, ".",
             f"ID={b.name}_{ref_allele}"]
        )
        # To show the label of all block in IGV
        # I break the relation (Remove parent attribute)
        if not igv_show_label:
            records[-1][-1] += f";Parent={gene_name}"
        pos += b.length

    # save
    with open(fname, "w") as f_gff:
        f_gff.write("##gff-version 3\n")
        f_gff.write("\n".join(["\t".join(record) for record in records]))
    return None


# load/save Genemsa
def load_msa(file_fasta: str, file_json: str) -> Genemsa:
    """
    load Genemsa from fasta and json

    Example:
      >>> from pyHLAMSA import Genemsa
      >>> a_gen.save_msa("a_gen.fa", "a_gen.json")
      >>> a_gen = Genemsa.load_msa("a_gen.fa", "a_gen.json")
    """
    # read
    msa = AlignIO.read(file_fasta, "fasta")
    with open(file_json) as f:
        msa.annotations = json.load(f)
    # main
    new_msa = Genemsa.from_MultipleSeqAlignment(msa)
    assert len(new_msa.get_reference()[1]) == new_msa.get_length()
    return new_msa


def save_msa(self: Genemsa, file_fasta: str, file_json: str):
    """ Save Genemsa to fasta and json """
    to_fasta(self, file_fasta, gap=True)
    with open(file_json, "w") as f:
        json.dump(self.meta_to_json(), f)
