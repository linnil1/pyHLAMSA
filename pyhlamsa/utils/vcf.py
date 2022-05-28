"""
VCF releated functions
"""
from datetime import datetime
from typing import List, Dict, Any
from collections import defaultdict
from dataclasses import dataclass, field
import pysam
from Bio import SeqIO

from . import cigar


@dataclass
class VcfVariant:
    """ The object to save the variant """
    pos: int
    ref: str
    alt: str
    chrom: str = ""

    def __post_init__(self):
        assert "-" not in self.ref
        assert "-" not in self.alt

    def __hash__(self):
        # Ignore chrom
        return hash((self.pos, self.ref, self.alt))


def extract_variants(ref_seq: str, tar_seq: str) -> List[VcfVariant]:
    """
    Compare two sequences and
    find the difference between them and
    format as VCF-like variants
    """
    pos = 0
    ref_gap = 0
    variants = []  # type: List[VcfVariant]
    prev_nomatch = False

    # remove when both sequences are gap
    bases = [i for i in zip(ref_seq, tar_seq) if not (i[0] == "-" and i[1] == "-")]
    ref_seq, tar_seq = map("".join, zip(*bases))

    # extract variant via cigar
    cigar_list = cigar.calculate_cigar(ref_seq, tar_seq)
    for cigar_op, cigar_count in cigar_list:
        if cigar_op == "M":
            pos += cigar_count
            prev_nomatch = False
        elif cigar_op == "X":
            if (prev_nomatch and len(variants[-1].ref) != len(variants[-1].alt)
                    and variants[-1].pos == 1):
                variants[-1] = VcfVariant(
                    pos=variants[-1].pos,
                    ref=variants[-1].ref,
                    alt=tar_seq[pos])
                pos += 1
                cigar_count -= 1
            for _ in range(cigar_count):
                variants.append(VcfVariant(pos=pos + 1 - ref_gap,
                                           ref=ref_seq[pos],
                                           alt=tar_seq[pos]))
                pos += 1
            prev_nomatch = True
        elif cigar_op == "I":
            if pos:
                if not prev_nomatch:
                    variants.append(VcfVariant(
                        pos=pos - ref_gap,
                        ref=ref_seq[pos - 1],
                        alt=tar_seq[pos - 1:pos + cigar_count]))
                else:
                    variants[-1] = VcfVariant(
                        pos=variants[-1].pos,
                        ref=variants[-1].ref,
                        alt=variants[-1].alt + tar_seq[pos:pos + cigar_count])
            else:  # pos == 0
                variants.append(VcfVariant(
                    pos=1,
                    ref=ref_seq[cigar_count],
                    alt=tar_seq[:cigar_count + 1]))
            pos += cigar_count
            ref_gap += cigar_count
            prev_nomatch = True
        elif cigar_op == "D":
            if pos:
                if not prev_nomatch:
                    variants.append(VcfVariant(
                        pos=pos - ref_gap,
                        ref=ref_seq[pos - 1:pos + cigar_count],
                        alt=tar_seq[pos - 1]))
                else:
                    variants[-1] = VcfVariant(
                        pos=variants[-1].pos,
                        ref=variants[-1].ref + ref_seq[pos:pos + cigar_count],
                        alt=variants[-1].alt)
            else:  # pos == 0
                variants.append(VcfVariant(
                    pos=1,
                    ref=ref_seq[:cigar_count + 1],
                    alt=tar_seq[cigar_count]))
            pos += cigar_count
            prev_nomatch = True
        else:
            raise ValueError(f"Invalid cigar {cigar_op}")

    assert len(ref_seq) == len(tar_seq) == pos
    return variants


def get_vcf_header(ref_name: str, ref_seq: str) -> str:
    """ Return minimal vcf header """
    return f"""\
##fileformat=VCFv4.1
##fileDate={datetime.now().strftime("%Y%m%d")}
##contig=<ID={ref_name},length={len(ref_seq.replace("-", ""))}>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"""


def variants_to_table(allele_variants: Dict[str, List[VcfVariant]]) -> List[List[Any]]:
    """ Summary variant in each allele """
    # Convert Dict[allele_name, List[variant]] = allele_variants
    # to      Dict[variant, List[allele_name]] = variant_alleles
    variant_alleles = defaultdict(set)
    for allele, variants in allele_variants.items():
        for variant in variants:
            variant_alleles[variant].add(allele)

    # collect all alleles_name in order
    alleles_name = sorted(allele_variants.keys())

    # set header
    table = []
    table.append(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO",
                  "FORMAT", *alleles_name])

    # each record
    for variant, alleles in variant_alleles.items():
        table.append([
            variant.chrom, variant.pos, ".", variant.ref, variant.alt,  # type: ignore
            ".", ".", ".",
            "GT", *("1" if name in alleles else "0" for name in alleles_name)
        ])
    return table


def read_vcf(file_vcf: str, file_fasta: str) -> Dict[str, Dict[str, str]]:
    """
    An experiment function to read vcf file into msa

    Only implemented:

    * haploid
    * single-allele
    * '.' not in GT
    * no-insertion
    """
    vcf = pysam.VariantFile(file_vcf)

    # create mutliple msa (dict type)
    # chrom_msa = {chrom_name: { allele_name: sequencs }}
    seqs = {seq.id: str(seq.seq) for seq in SeqIO.parse(file_fasta, "fasta")}
    chrom_msa = {str(name): {str(name): seqs[name]} for name in vcf.header.contigs}
    for sample_name in vcf.header.samples:
        for ref, seq_dict in chrom_msa.items():
            seq_dict[str(sample_name)] = seq_dict[ref]

    # modify msa by variant
    for record in vcf.fetch():
        # not implement mutliple alts with GT
        assert record.ref
        assert record.alts
        assert len(record.alts) == 1
        ref = record.ref
        alt = record.alts[0]
        # not implement insertion
        assert len(ref) >= len(alt)

        for allele_name, sample in record.samples.items():
            allele_name = str(allele_name)
            # only implement haploid
            assert "GT" in sample
            assert len(sample["GT"]) == 1
            gt = str(sample["GT"][0])
            assert gt in "01"
            if gt == "1":
                sample_seq = chrom_msa[record.chrom][allele_name]
                pos = record.pos - 1
                # must be True if vcf is correct
                assert sample_seq[pos:pos + len(ref)] == ref
                if pos:
                    alt_with_gap = alt + "-" * (len(ref) - len(alt))
                else:
                    alt_with_gap = "-" * (len(ref) - len(alt)) + alt
                sample_seq = \
                    sample_seq[:pos] \
                    + alt_with_gap \
                    + sample_seq[pos + len(ref):]
                chrom_msa[record.chrom][allele_name] = sample_seq

    return chrom_msa
