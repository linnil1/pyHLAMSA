import logging
import os
from glob import glob
import re
from pprint import pprint

logging.basicConfig(level=logging.INFO)


class HLAmsa:
    def __init__(self, genes=[], filetype=["gen", "nuc"],
                 imgt_folder="alignments", version="3430"):
        self.imgt_folder = imgt_folder
        if not os.path.exists(self.imgt_folder):
            self.download(True, version=version)
        else:
            logging.info(f"IMGT ver={version} exists")
        assert os.path.exists(self.imgt_folder)

        if not genes:
            fs = glob(f"{self.imgt_folder}/*.txt")
            genes = set([f.split("/")[-1].split("_")[0] for f in fs])

        self.genes = {}
        for gene in genes:
            logging.info(f"Reading {gene}")
            self.genes[gene] = self.parse_alignment(f"{self.imgt_folder}/{gene}_gen.txt")

    def download(self, download=True, version="3430"):
        """
        Download the IMGTHLA alignments folder to `imgt_folder`
        """
        if os.path.exists(self.imgt_folder):
            return
        # TODO: Auto find the latest version
        fname = f"Alignments_Rel_{version}"
        url = "ftp://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/" + fname + ".zip"
        logging.info(f"Download IMGT data from {url}")
        os.system(f"wget {url}")
        os.system(f"unzip {fname}.zip")
        os.system(f"mv alignments {self.imgt_folder}")

    def parse_alignment(self, fname):
        """
        Read `fname` alignments files with IMGT-alignment format
        """

        # parse aligments
        alleles = {}
        ref_allele = ""
        for line in open(fname):
            line = line.strip()
            if not re.match(r"^\w+\*", line):
                continue

            match = re.findall(r"^(.*?) +(.*)", line)
            assert match
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
            alleles[allele] = ''.join(seq)

        for allele in rm_allele:
            logging.warning(f"Remove {allele} due to length in {fname}")
            del alleles[allele]
        return alleles


if __name__ == "__main__":
    msa = HLAmsa(["A"], imgt_folder="alignments")

