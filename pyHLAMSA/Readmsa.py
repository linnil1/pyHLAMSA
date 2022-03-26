import re
import logging
from typing import List, Tuple, Dict, Tuple
from Bio import AlignIO, SeqIO
from Bio.Align import MultipleSeqAlignment
from . import Genemsa

logger = logging.getLogger(__name__)


def from_MSF_file(file_msf: str) -> Genemsa:
    """ Read .msf file """
    return Genemsa.from_MultipleSeqAlignment(AlignIO.read(file_msf, "msf"))


def from_alignment_file(fname: str, seq_type="") -> Genemsa:
    """ Read MSA format file `fname` .e.g. IMGT/alignments/A_gen.txt """
    # read all alleles
    alleles = read_alignment_sequence(fname)
    new_msa = Genemsa("")
    new_msa.alleles = {name: seq.replace("|", "") for name, seq in alleles.items()}

    # use first sequence as reference
    ref_seq = list(alleles.values())[0]
    new_msa.blocks = [{'length': len(seq)} for seq in ref_seq.split("|")]

    if seq_type:
        new_msa.seq_type = seq_type
        new_msa._assume_label()

    # check
    leng = new_msa.get_length()
    for name in new_msa.get_sequence_names():
        assert len(new_msa.get(name)) == leng
    return new_msa


def read_alignment_sequence(fname: str) -> Dict[str, str]:
    """
    Read MSA format file `fname`

    Returns:
        alleles (dict of str,str): The dictionary that map allele name to sequence
    """
    # parse aligments
    alleles = {}
    ref_allele = ""
    for line in open(fname):
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
        logger.warning(f"Remove {allele} due to length in {fname}")
        del alleles[allele]
    return alleles


def read_dat_block(file_dat: str):
    """
    Read block information in .dat file

    Args:
      file_dat (str): File path to xx.dat

    Returns:
      dat_object(dict): {allele_name: [[block_name(str), start(int), end(int)],]
          , where block_name is intron1, exon1, ...
          , and the start, end is the position without gap
    """
    now_allele = ""
    read_next = False
    data = {}

    for line in open(file_dat):
        # read allele name
        if "allele=" in line:
            now_allele = line.split('"')[1]
            now_allele = now_allele.replace("HLA-", "")
            # no duplicated allele_name
            assert now_allele not in data
            data[now_allele] = []

        # read blocks
        elif re.match(r"FT\s+((UTR)|(exon)|(intron))", line):
            line = line.split()
            data[now_allele].append([line[1], *map(lambda a: int(a), line[2].split(".."))])
            if line[1] != "UTR":
                read_next = True

        # read block name .e.g. exon1 exon2
        elif read_next:
            read_next = False
            data[now_allele][-1][0] += line.split('"')[1]

    # make sure the order start and end position in dat are correct
    for allele in data:
        data[allele] = list(sorted(data[allele], key=lambda a: a[1]))

    return data


def apply_dat_info_on_msa(msa, dat) -> Genemsa:
    """
    Merge dat object

    Args:
        dat (object): Object created from `Genemsa.read_dat`

    Returns:
        Genemsa
    """
    if not len(msa):
        raise ValueError("MSA is empty")
    msf_length = len(msa._get_first()[1])
    if msa.seq_type != "gen" and msa.seq_type != "nuc":
        raise TypeError("Check this object is gen or nuc")
    dat = copy.deepcopy(dat)
    # I don't know how it work...
    # need more documentation

    # Check consistent between dat and msf
    block_cord = {}
    new_alleles = {}
    for allele, seq in msa.alleles.items():
        # Check allele name in the dat
        # and sequence is consistent to dat
        hla_name = allele
        assert hla_name in dat

        # extract exon region
        if msa.seq_type == "nuc":
            dat_hla = []
            end = 0
            for i in [i for i in dat[hla_name] if "exon" in i[0]]:
                gap = i[1] - end - 1
                dat_hla.append([i[0], i[1] - gap, i[2] - gap])
                end = dat_hla[-1][2]
            dat[hla_name] = dat_hla

        if len(seq.replace("-", "")) != dat[hla_name][-1][2]:
            logger.warning(f"Ignore {allele}, msf length is different from dat")
            continue

        # rename UTR
        # assume blocks are ordered
        if "UTR" == dat[hla_name][0][0]:
            dat[hla_name][0][0] = "5UTR"
        if "UTR" == dat[hla_name][-1][0]:
            dat[hla_name][-1][0] = "3UTR"

        # Coordination mapping from sequence to msf
        seq_cord = []
        assert len(seq) == msf_length
        for i in range(len(seq)):
            if seq[i] != "-":
                assert seq[i] in "ATCG"
                seq_cord.append(i)

        # Get the position of each block in msf
        old_cord = copy.deepcopy(block_cord)
        for block_name, start, end in dat[hla_name]:
            if block_name not in block_cord:
                block_cord[block_name] = [99999, 0]

            block_cord[block_name][0] = min(block_cord[block_name][0],
                                            seq_cord[start - 1])
            block_cord[block_name][1] = max(block_cord[block_name][1],
                                            seq_cord[end - 1] + 1)

        # Check consistent
        # Assume first allele has correct coordination
        # , otherwise all allele will ignore
        block_cord_list = list(sorted(block_cord.values()))
        for i in range(1, len(block_cord_list)):
            if block_cord_list[i][0] < block_cord_list[i - 1][1]:
                logger.warning(f"Ignore {allele}, msf inconsistent")
                block_cord = old_cord
                break
        else:
            new_alleles[allele] = seq

    # Reset position of each block
    # If perfect -> start_i == end_{i-1}
    # If not -> start_i = end_{i-1}
    block_cord_list = list(sorted(block_cord.items(), key=lambda i: i[1]))
    # TODO: is end == msf_length
    end = msf_length
    # end = block_cord_list[-1][1][1]
    for i in reversed(range(len(block_cord_list))):
        block_cord_list[i], end = ([block_cord_list[i][0],
                                    [block_cord_list[i][1][0], end]],
                                   block_cord_list[i][1][0])

    # Check the sequence length in msf is same as block coordination
    alleles = new_alleles
    new_alleles = {}
    block_cord = dict(block_cord_list)
    for allele, seq in alleles.items():
        allele_block = {i[0]: i[2] - i[1] + 1 for i in dat[allele]}

        for block_name in block_cord:
            len_block = len(seq[block_cord[block_name][0]:block_cord[block_name][1]].replace("-", ""))
            len_block_ref = allele_block.get(block_name, 0)
            # should always true, otherwise ignore in above
            assert len_block == len_block_ref
        else:
            new_alleles[allele] = seq

    # To new object
    new_msa = msa.copy()
    new_msa.alleles = new_alleles
    new_msa.blocks = []
    for name, (start, end) in block_cord_list:
        b = {'length': end - start}
        if "3UTR" == name:
            b.update({'type': "three_prime_UTR", 'name': "3UTR"})
        elif "5UTR" == name:
            b.update({'type': "five_prime_UTR", 'name': "5UTR"})
        elif "exon" in name:
            b.update({'type': "exon", 'name': name})
        elif "intron" in name:
            b.update({'type': "intron", 'name': name})
        new_msa.blocks.append(b)
    return new_msa
