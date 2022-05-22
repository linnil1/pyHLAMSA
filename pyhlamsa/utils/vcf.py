"""
VCF releated functions
"""
from datetime import datetime
from typing import List, Dict, Any
from collections import defaultdict
from dataclasses import dataclass, field

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
    variants = []
    prev_indel = False

    # remove when both sequences are gap
    bases = [i for i in zip(ref_seq, tar_seq) if not (i[0] == "-" and i[1] == "-")]
    ref_seq, tar_seq = map(''.join, zip(*bases))

    # extract variant via cigar
    cigar_list = cigar.calculate_cigar(ref_seq, tar_seq)
    for cigar_op, cigar_count in cigar_list:
        if cigar_op == "M":
            pos += cigar_count
            prev_indel = False
        elif cigar_op == "X":
            for i in range(cigar_count):
                variants.append(VcfVariant(pos=pos + 1 - ref_gap,
                                           ref=ref_seq[pos],
                                           alt=tar_seq[pos]))
                pos += 1
            prev_indel = False
        elif cigar_op == "I":
            if pos:
                if not prev_indel:
                    variants.append(VcfVariant(
                        pos=pos - ref_gap,
                        ref=ref_seq[pos - 1],
                        alt=tar_seq[pos - 1:pos + cigar_count]))
                else:
                    variants.append(VcfVariant(
                        pos=variants[-1].pos,
                        ref=variants[-1].ref,
                        alt=variants[-1].alt + tar_seq[pos:pos + cigar_count]))
            else:  # pos == 0
                variants.append(VcfVariant(
                    pos=1,
                    ref=ref_seq[cigar_count],
                    alt=tar_seq[:cigar_count + 1]))
            pos += cigar_count
            ref_gap += cigar_count
            prev_indel = True
        elif cigar_op == "D":
            if pos:
                if not prev_indel:
                    variants.append(VcfVariant(
                        pos=pos - ref_gap,
                        ref=ref_seq[pos - 1:pos + cigar_count],
                        alt=tar_seq[pos - 1]))
                else:
                    variants.append(VcfVariant(
                        pos=variants[-1].pos,
                        ref=variants[-1].ref + ref_seq[pos:pos + cigar_count],
                        alt=variants[-1].alt))
            else:  # pos == 0
                variants.append(VcfVariant(
                    pos=1,
                    ref=ref_seq[:cigar_count + 1],
                    alt=tar_seq[cigar_count]))
            pos += cigar_count
            prev_indel = True
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
            variant.chrom, variant.pos, '.', variant.ref, variant.alt,  # type: ignore
            '.', '.', '.',
            'GT', *('1' if name in alleles else '0' for name in alleles_name)
        ])
    return table
