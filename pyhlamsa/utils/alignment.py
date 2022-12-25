"""IMGT-alignment format (xx_gen.txt) related code"""
import re
import logging

logger = logging.getLogger(__name__)


def parse_alignment_txt(fname: str) -> dict[str, str]:
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
