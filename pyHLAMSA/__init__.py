from __future__ import annotations
import logging
import os
from glob import glob
from pprint import pprint
from typing import List, Tuple, Dict

from .Genemsa import Genemsa


class HLAmsa:
    """
    A HLA interface

    Attributes:
        genes (dict of str, Genemsa): The msa object for each gene
    """

    def __init__(self, genes=[], filetype=["gen", "nuc"],
                 imgt_folder="alignments", version="3430"):
        """
        The instance will download the IMGT alignment folder into `imgt_folder`
        with version `verion` and read the `{gene}_{filetype}.txt` files.

        Args:
            genes (list of str): A list of genes you want to read.

                Leave Empty if you want read all gene in HLA

                set None to read none of gene.

            filetype (list of str): A list of filetype.

                If both `gen` and `nuc` are given, it will merge them automatically.
        """
        self.logger = logging.getLogger(__name__)
        self.imgt_folder = imgt_folder
        if not os.path.exists(self.imgt_folder):
            self._download(True, version=version)
        else:
            self.logger.info(f"IMGT ver={version} exists")
        assert os.path.exists(self.imgt_folder)

        # gene
        self.genes = {}
        if genes is None:
            return
        if not genes:
            genes = self.list_gene()

        self.logger.info(f"Read Gene {genes}")
        for gene in genes:
            self.logger.info(f"Reading {gene}")
            self.genes[gene] = self._read_alignments(gene, filetype)
            self.logger.debug(f"Merged {self.genes[gene]}")

    def list_gene(self) -> List[str]:
        """ List the gene in folder """
        fs = glob(f"{self.imgt_folder}/*.txt")
        genes = set([f.split("/")[-1].split("_")[0] for f in fs])
        # TODO:
        # * P -> only exon
        # * E last exon not shown in nuc
        # * N the nuc has one more bp than in gen
        genes = genes - set(["ClassI", "DRB", "N", "P", "E"])
        genes = sorted(list(genes))
        return genes

    def _read_alignments(self, gene: str, filetype: List[str]):
        """
        Read `{gene}_{filetype}.txt`.

        If both `gen` and `nuc` are given, it will merge them.

        The alignment data is stored in `self.genes[gene]` in `Genemsa` instance
        """
        gene_merge_exclude = []
        if "gen" in filetype:
            msa_gen = Genemsa(gene, seq_type="gen")
            msa_gen.read_alignment_file(f"{self.imgt_folder}/{gene}_gen.txt")
            self.logger.debug(f"{msa_gen}")
        if "nuc" in filetype and gene not in gene_merge_exclude:
            msa_nuc = Genemsa(gene, seq_type="nuc")
            # Special Case: DRB* nuc are in DRB_nuc.txt
            if gene.startswith("DRB"):
                msa_nuc.read_alignment_file(f"{self.imgt_folder}/DRB_nuc.txt")
                self.logger.debug(f"DRB: {msa_nuc}")
                msa_nuc = msa_nuc.select_allele(gene + ".*")
                self.logger.debug(f"{msa_nuc}")
            else:
                msa_nuc.read_alignment_file(f"{self.imgt_folder}/{gene}_nuc.txt")
                self.logger.debug(f"{msa_nuc}")

        if "gen" in filetype and "nuc" in filetype \
                and gene not in gene_merge_exclude:
            return msa_gen.merge_exon(msa_nuc)
        elif "gen" in filetype:
            return msa_gen
        elif "nuc" in filetype:
            return msa_nuc
        else:
            return None

    def _download(self, download=True, version="3430"):
        """
        Download the IMGTHLA alignments folder to `imgt_folder`
        """
        if os.path.exists(self.imgt_folder):
            return
        # TODO: Auto find the latest version
        fname = f"Alignments_Rel_{version}"
        url = "ftp://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/" + fname + ".zip"
        self.logger.info(f"Download IMGT data from {url}")
        os.system(f"wget {url}")
        os.system(f"unzip {fname}.zip")
        os.system(f"mv alignments {self.imgt_folder}")


class HLAmsaEX(HLAmsa):
    """
    A HLA interface but read gene from MSF and hla.dat

    Attributes:
        genes (dict of str, Genemsa): The msa object for each gene
    """

    def __init__(self, genes=[], filetype=["gen", "nuc"],
                 imgt_folder="IMGT", version="3430"):
        """
        The instance will download the IMGT alignment folder into `imgt_folder`
        with version `verion` and read the `{gene}_{filetype}.msf` files.

        Args:
            genes (list of str): A list of genes you want to read.

                Leave Empty if you want read all gene in HLA

                set None to read none of gene.

            filetype (list of str): A list of filetype.

                If both `gen` and `nuc` are given, it will merge them automatically.
        """
        self.logger = logging.getLogger(__name__)
        self.imgt_folder = imgt_folder
        if not os.path.exists(self.imgt_folder):
            self._download(True, version=version)
        else:
            self.logger.info(f"IMGT ver={version} exists")
        assert os.path.exists(self.imgt_folder)

        # gene
        self.genes = {}
        if genes is None:
            return
        if not genes:
            genes = self.list_gene()

        self.logger.info(f"Read Gene {genes}")
        for gene in genes:
            self.logger.info(f"Reading {gene}")
            self.genes[gene] = self._read_msf(gene, filetype)
            self.logger.debug(f"Merged {self.genes[gene]}")

    def list_gene(self) -> List[str]:
        """ List the gene in folder """
        fs = glob(f"{self.imgt_folder}/msf/*.msf")
        genes = set([f.split("/")[-1].split("_")[0] for f in fs])
        # E_nuc, S_nuc, N_nuc: exon is different from nuc
        # MICA MICB TAP1 TAP2 has no record in hla.dat
        genes = genes - set(["DRB", "DRB345", "MICA", "MICB", "TAP1", "TAP2",
                             "E", "S", "N"])
        genes = sorted(list(genes))
        return genes

    def _read_msf(self, gene: str, filetype: List[str]):
        """
        Read `{gene}_{filetype}.txt`.

        If both `gen` and `nuc` are given, it will merge them.

        The alignment data is stored in `self.genes[gene]` in `Genemsa` instance
        """
        self.dat = Genemsa.read_dat(f"{self.imgt_folder}/hla.dat")

        gene_merge_exclude = ["P"]
        if "gen" in filetype:
            msa_gen = Genemsa(gene, seq_type="gen")
            msa_gen.read_MSF_file(f"{self.imgt_folder}/msf/{gene}_gen.msf")
            msa_gen = msa_gen.merge_dat(self.dat)
            self.logger.debug(f"{msa_gen}")

        if "nuc" in filetype and gene not in gene_merge_exclude:
            msa_nuc = Genemsa(gene, seq_type="nuc")
            # Special Case: DRB* nuc are in DRB_nuc.txt
            if gene.startswith("DRB"):
                msa_nuc.read_MSF_file(f"{self.imgt_folder}/msf/DRB_nuc.msf")
                msa_nuc = msa_nuc.select_allele(gene + ".*")
            else:
                msa_nuc.read_MSF_file(f"{self.imgt_folder}/msf/{gene}_nuc.msf")
            msa_nuc = msa_nuc.merge_dat(self.dat)
            self.logger.debug(f"{msa_nuc}")

        if "gen" in filetype and "nuc" in filetype \
                and gene not in gene_merge_exclude:
            return msa_gen.merge_exon(msa_nuc)
        elif "gen" in filetype:
            return msa_gen
        elif "nuc" in filetype:
            return msa_nuc
        else:
            return None

    def _download(self, download=True, version="3430"):
        """
        Download the IMGTHLA alignments folder to `imgt_folder`
        """
        if os.path.exists(self.imgt_folder):
            return
        # TODO: Auto find the latest version
        os.system(f"git clone https://github.com/ANHIG/IMGTHLA.git {self.imgt_folder}")
        os.system(f"cd {self.imgt_folder} && git checkout {version} && cd ..")
        os.system(f"cd {self.imgt_folder} && git lfs pull && cd ..")
