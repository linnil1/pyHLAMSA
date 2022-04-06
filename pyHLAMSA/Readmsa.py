import re
import copy
import logging
from collections import defaultdict
from typing import List, Tuple, Dict, Tuple
from Bio import AlignIO, SeqIO
from Bio.Align import MultipleSeqAlignment
from . import Genemsa, BlockInfo

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
    new_msa.blocks = [BlockInfo(length=len(seq)) for seq in ref_seq.split("|")]

    if seq_type:
        new_msa.seq_type = seq_type
        new_msa._assume_label()  # include label

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
        logger.warning(f"{allele} violate msa length in {fname}. Removed")
        del alleles[allele]
    return alleles


def read_dat_block(file_dat: str) -> Dict:
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
        if "/allele=" in line:
            now_allele = line.split('"')[1]
            now_allele = now_allele.replace("HLA-", "")
            # no duplicated allele_name
            assert now_allele not in data
            data[now_allele] = []
        elif not line.startswith("FT"):
            continue

        # read blocks
        elif re.match(r"FT\s+((UTR)|(exon)|(intron))", line):
            line = line.split()
            start, end = list(map(lambda a: int(a), line[2].split("..")))
            data[now_allele].append({
                'name': line[1],
                'start': start,
                'end': end,
            })

        # read block name .e.g. exon1 exon2
        elif "/number=" in line:
            data[now_allele][-1]['name'] += line.split('"')[1]
        elif "/pseudo" in line:
            data[now_allele][-1]['pseudo'] = True
        else:
            continue

    # make sure the order start and end position in dat are correct
    for allele in data:
        data[allele] = list(sorted(data[allele], key=lambda a: a['start']))
    return data


def apply_dat_info_on_msa(msa: Genemsa, dat: Dict) -> Genemsa:
    """
    Apply the dat information to MSA to cut the exon, intron position.

    This function will estimate the intron/exon position in MSA,
    when given the intron/exon position in no-gap sequence

    Args:
        dat (object): Object created from `Genemsa.read_dat`

    Returns:
        Genemsa
    """
    if not len(msa):
        raise ValueError("MSA is empty")
    msf_length = len(msa.get_reference()[1])
    if msa.seq_type != "gen" and msa.seq_type != "nuc":
        raise TypeError("Check this object is gen or nuc")
    dat = copy.deepcopy(dat)

    # block_cord save the min and max possible position of intron and exon
    # first element: the maximum of the first position
    # first element: the minimum of the last position
    block_cord = defaultdict(lambda: [msf_length, 0])
    # new_alleles: save the sequence
    new_alleles = {}

    for allele_name, seq in msa.alleles.items():
        if allele_name not in dat:
            logger.warning(f"Ignore {allele_name}: not exist in dat")
            continue

        # make sure the intron/exon order
        dat[allele_name] = sorted(dat[allele_name], key=lambda i: i['start'])

        # rename UTR
        if "UTR" == dat[allele_name][0]['name']:
            dat[allele_name][0]['name'] = "5UTR"
        if "UTR" == dat[allele_name][-1]['name']:
            dat[allele_name][-1]['name'] = "3UTR"

        # type is nuc i.e. no exon exist
        # remove all intron, and recalculate the position
        if msa.seq_type == "nuc":
            dat_hla = []
            end = 0
            for d in filter(lambda i: "exon" in i['name'], dat[allele_name]):
                # set pseduo gene length = 0
                if d.get('pseudo', False):
                    d['end'] = d['start'] - 1
                # gap  = total previous intron length
                gap = d['start'] - end - 1
                dat_hla.append({
                    'name': d['name'],
                    'start': d['start'] - gap,
                    'end': d['end'] - gap
                })
                end = dat_hla[-1]['end']
            dat[allele_name] = dat_hla

        # check sequence length
        if len(seq.replace("-", "")) != dat[allele_name][-1]['end']:
            logger.warning(f"Ignore {allele_name}, "
                           f"msf length is different from dat. "
                           f"seq: {len(seq.replace('-', ''))} "
                           f"dat: {dat[allele_name][-1]['end']}")
            continue
        # assert len(seq) == msf_length
        # check sequence content
        if not all(i in "ATCG-" for i in seq):
            logger.warning(f"Ignore {allele_name}, msf length contains invalid character"
                           " i.e. character shoudle be one of ATCG-(gap)")
            continue

        # Coordination mapping from sequence to gapped-sequence
        seq_cord = []
        for i in range(len(seq)):
            if seq[i] != "-":
                seq_cord.append(i)

        # Calculate the intron/exon position after insert this allele
        # If the intron/exon overlap -> ignore this allele, restore the position
        old_cord = copy.deepcopy(block_cord)
        for d in dat[allele_name]:
            block_name, start, end = d['name'], d['start'], d['end']
            # calculate the position
            block_cord[block_name][0] = min(block_cord[block_name][0],
                                            seq_cord[start - 1])
            block_cord[block_name][1] = max(block_cord[block_name][1],
                                            seq_cord[end - 1] + 1)
        block_cord_list = list(sorted(block_cord.values()))
        for i in range(1, len(block_cord_list)):
            # is overlap
            if block_cord_list[i - 1][1] > block_cord_list[i][0]:
                logger.warning(f"Ignore {allele_name}, maybe something wrong in dat")
                block_cord = old_cord
                break
        else:
            # success
            new_alleles[allele_name] = seq

    # Because the ensure two consecutive region will not overlap before,
    # we just set any value between the lower and upper bound
    # of each intron/exon position in block_cord
    block_cord_list = []
    start_pos = 0
    for block_name, (start, end) in sorted(block_cord.items(), key=lambda i: i[1]):
        assert start_pos <= end
        block_cord_list.append((block_name, [start_pos, end]))
        start_pos = end
    assert block_cord_list[-1][1][1] <= msf_length
    block_cord_list[-1][1][1] = msf_length

    # double check the regions is match dat information
    alleles = new_alleles
    new_alleles = {}
    for allele_name, seq in alleles.items():
        block_length = {i['name']: i['end'] - i['start'] + 1 for i in dat[allele_name]}

        for block_name, (start, end) in block_cord_list:
            len_block = len(seq[start:end].replace("-", ""))
            len_block_ref = block_length.get(block_name, 0)
            assert len_block == len_block_ref
        else:
            new_alleles[allele_name] = seq

    # New msa
    new_msa = msa.copy()
    new_msa.blocks = []
    for name, (start, end) in block_cord_list:
        new_msa.blocks.append(BlockInfo(
            length=end - start,
            name=name,
            type={
                "3UTR": "three_prime_UTR",
                "5UTR": "five_prime_UTR",
                "exon": "exon",
                "intron": "intron",
            }.get(name),
        ))
    new_msa = new_msa.reset_index()
    new_msa.alleles = new_alleles
    return new_msa
