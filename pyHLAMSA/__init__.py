from __future__ import annotations
import os
import logging
import subprocess
from glob import glob
from typing import List, Set

from .Genemsa import Genemsa
from . import Readmsa


class Familymsa:
    """
    A abstract class to handle Gene Family
    """
    def __init__(self, genes=[], filetype=["gen", "nuc"],
                 db_folder="alignments", version="latest"):
        """
        When init, it will

        * Init parameter
        * Download db
        * Read db to Genemsa

        Args:
            genes (str or list of str): A list of genes you want to read.

                Leave Empty if you want read all gene in HLA

                set None to read none of gene.

            filetype (str or list of str): A list of filetype.

                If both `gen` and `nuc` are given, it will merge them automatically.
        """
        self.logger = logging.getLogger(__name__)
        self.db_folder = db_folder
        self.genes = {}

        if genes is None:
            return
        if type(genes) is str:
            genes = [genes]
        if type(filetype) is str:
            filetype = [filetype]

        filetype = set(filetype)
        if not genes:  # if empty -> read all
            genes = self._list_db_gene(filetype)

        # main
        self._download(version)
        for gene_name in genes:
            self.logger.info(f"Reading {gene_name}")
            self.genes[gene_name] = self._read_db_gene(gene_name, filetype)

    def list_genes(self):
        """ List all the gene's name in this family """
        return list(self.genes.keys())

    def __getitem__(self, index: str) -> Genemsa:
        """ Get specific gene's msa """
        return self.genes[index]

    def _download(self, version: str):
        """ Check before running download """
        if not os.path.exists(self.db_folder):
            self._download_db(version=version)
            if not os.path.exists(self.db_folder):
                raise ValueError("Fail to download")
        else:
            self.logger.info(f"{self.db_folder} exists")

    def _download_db(self, version):
        """ Abstract method: code for downloading your db """
        raise NotImplementedError

    def _run_shell(self, *args, cwd=None):
        """ Run shell code """
        self.logger.debug("Run " + " ".join(args))
        with subprocess.Popen(args, cwd=cwd) as proc:
            proc.wait()

    def _list_db_gene(self, filetype=set(["gen", "nuc"])) -> List[str]:
        """ Abstract method: code for listing gene names """
        raise NotImplementedError

    def _read_db_gene(self, gene, filetype=set(["gen", "nuc"])):
        """ Abstract method: code for reading and merging function """
        raise NotImplementedError


class HLAmsa(Familymsa):
    """
    A HLA interface

    Attributes:
        genes (dict of str, Genemsa): The msa object for each gene
    """

    def __init__(self, genes=[], filetype=["gen", "nuc"],
                 imgt_alignment_folder="alignments_v3430", version="3430"):
        super().__init__(genes, filetype,
                         db_folder=imgt_alignment_folder, version=version)

    def _download_db(self, version: str):
        """ Download the IMGTHLA alignments folder to `db_folder` """
        # TODO: Auto find the latest version
        fname = f"Alignments_Rel_{version}"
        url = "ftp://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/" + fname + ".zip"
        self.logger.info(f"Download IMGT data from {url}")
        self._run_shell("wget", url)
        self._run_shell("unzip", f"{fname}.zip")
        self._run_shell("mv", "alignments", self.db_folder)

    def _get_name(self, search_name: str) -> Set[str]:
        """ Handy function to list names from file pattern """
        arr_files = glob(search_name)
        return set([f.split("/")[-1].split("_")[0] for f in arr_files])

    def _list_db_gene(self, filetype=set(["gen", "nuc"])) -> List[str]:
        """ List the gene in folder """
        DRB = set(['DRB1', 'DRB3', 'DRB4', 'DRB5'])
        if "gen" in filetype:
            names = names_gen = self._get_name(f"{self.db_folder}/*_gen.txt") | DRB
        if "nuc" in filetype:
            names = names_nuc = self._get_name(f"{self.db_folder}/*_nuc.txt") | DRB
        if "gen" in filetype and "nuc" in filetype:
            names = names_gen & names_nuc
            # * P -> only exon
            # * N the nuc has one more bp than in gen
            #    <N gen alleles=5 block=5UTR(225) exon1(161) 3UTR(250)>
            #    <N nuc alleles=5 block=exon1(162)>
            # * E exon7 is ambiguous, exon8 is gone
            #    <E gen alleles=258 block=5UTR(301) exon1(64) intron1(130) exon2(270) intron2(244) exon3(276) intron3(621) exon4(276) intron4(124) exon5(117) intron5(751) exon6(33) intron6(104) exon7(43) intron7(165) exon8(5) 3UTR(356)>
            #    <E nuc alleles=262 block=exon1(64) exon2(270) exon3(276) exon4(276) exon5(117) exon6(33) exon7(41)>
            # * S exon1 is ambiguous
            #   <S gen alleles=0 block=5UTR(222) exon1(28) intron1(104) exon2(48) intron2(161) exon3(192) 3UTR(120)>
            #   <S gen alleles=7 block=5UTR(222) exon1(27) intron1(104) exon2(48) intron2(161) exon3(192) 3UTR(120)>
            if "gen" in filetype and "nuc" in filetype:
                names = names - set(["N", "E", "S"])
        return list(sorted(names))

    def _read_db_gene(self, gene: str, filetype=set(["gen", "nuc"])):
        """
        Read `{gene}_{filetype}.txt`.

        If both `gen` and `nuc` are given, it will merge them.
        """
        if "gen" in filetype:
            msa_gen = Readmsa.from_alignment_file(f"{self.db_folder}/{gene}_gen.txt")
            msa_gen.seq_type = "gen"
            msa_gen.gene_name = gene
            if gene != "P":
                # P: special case: has even block
                # <P gen alleles=5 block=(475) (261) (589) (276) (124) (117) (412) (33) (150) (48) (163) (297)>
                msa_gen._assume_label()
            self.logger.debug(f"Read {msa_gen}")
        if "nuc" in filetype:
            # Special Case: DRB* nuc are in DRB_nuc.txt
            if gene.startswith("DRB"):
                msa_nuc = Readmsa.from_alignment_file(f"{self.db_folder}/DRB_nuc.txt")
                msa_nuc = msa_nuc.select_allele(gene + ".*")
                self.logger.debug(f"{msa_nuc}")
            else:
                msa_nuc = Readmsa.from_alignment_file(f"{self.db_folder}/{gene}_nuc.txt")
                self.logger.debug(f"{msa_nuc}")

            msa_nuc.seq_type = "nuc"
            msa_nuc.gene_name = gene
            msa_nuc._assume_label()

        if "gen" in filetype and "nuc" in filetype:
            # remove some gene
            diff_name = list(set(msa_gen.get_sequence_names()) - set(msa_nuc.get_sequence_names()))
            if diff_name:
                self.logger.warning(f"Remove alleles doesn't exist in gen and nuc either: {diff_name}")
            msa_gen = msa_gen.remove(diff_name)

            # merge
            msa_merged = msa_gen.merge_exon(msa_nuc)
            self.logger.debug(f"{msa_merged}")
            return msa_merged

        elif "gen" in filetype:
            return msa_gen
        elif "nuc" in filetype:
            return msa_nuc
        else:
            return None


class HLAmsaEX(Familymsa):
    """
    A HLA interface but read gene from MSF and
    read intron/exon information from hla.dat

    I think this one is much reliable

    Attributes:
        genes (dict of str, Genemsa): The msa object for each gene
    """

    def __init__(self, genes=[], filetype=["gen", "nuc"],
                 imgt_folder="IMGT_v3430", version="3430"):
        """
        The instance will download the IMGT/HLA into `imgt_folder`
        with version `verion` and read the `msf/{gene}_{filetype}.msf` files.

        Args:
            genes (list of str): A list of genes you want to read.

                Leave Empty if you want read all gene in HLA

                set None to read none of gene.

            filetype (list of str): A list of filetype.

                If both `gen` and `nuc` are given, it will merge them automatically.
        """
        super().__init__(genes, filetype,
                         db_folder=imgt_folder, version=version)

    def _download_db(self, download=True, version="3430"):
        """
        Download the IMGT/HLA msf and hla.dat to folder `imgt_folder`
        """
        # TODO: Auto find the latest version
        self._run_shell("git", "clone", "https://github.com/ANHIG/IMGTHLA.git", self.db_folder)
        self._run_shell("git", "checkout", version, cwd=self.db_folder)
        self._run_shell("git", "lfs", "pull", cwd=self.db_folder)

    def _get_name(self, search_name: str):
        """ Extract name from file pattern """
        arr_files = glob(search_name)
        return set([f.split("/")[-1].split("_")[0] for f in arr_files])

    def _list_db_gene(self, filetype) -> List[str]:
        """ List the gene in folder """
        DRB = set(['DRB1', 'DRB3', 'DRB4', 'DRB5'])
        if "gen" in filetype:
            names = names_gen = self._get_name(f"{self.db_folder}/msf/*_gen.msf") | DRB
        if "nuc" in filetype:
            # Most's of E nuc is different from dat record
            names = names_nuc = (self._get_name(f"{self.db_folder}/msf/*_nuc.msf") | DRB) - set(["E"])
        if "gen" in filetype and "nuc" in filetype:
            names = names_gen & names_nuc
        return sorted(names)

    def _read_db_gene(self, gene: str, filetype: List[str]):
        """
        Read `msf/{gene}_{filetype}.msf`.

        If both `gen` and `nuc` are given, it will merge them.
        """
        if not hasattr(self, "dat"):
            self.logger.debug(f"Reading hla.dat")
            self.dat = Readmsa.read_dat_block(f"{self.db_folder}/hla.dat")

        if "gen" in filetype:
            msa_gen = Readmsa.from_MSF_file(f"{self.db_folder}/msf/{gene}_gen.msf")
            msa_gen.seq_type = "gen"
            msa_gen.gene_name = gene
            msa_gen = Readmsa.apply_dat_info_on_msa(msa_gen, self.dat)
            self.logger.debug(f"{msa_gen}")

        if "nuc" in filetype:
            msa_nuc = Genemsa(gene, seq_type="nuc")
            # Special Case: DRB* nuc are in DRB_nuc.txt
            if gene.startswith("DRB"):
                msa_nuc = Readmsa.from_MSF_file(f"{self.db_folder}/msf/DRB_nuc.msf")
                msa_nuc = msa_nuc.select_allele(gene + ".*")
            else:
                msa_nuc = Readmsa.from_MSF_file(f"{self.db_folder}/msf/{gene}_nuc.msf")
            msa_nuc.seq_type = "nuc"
            msa_nuc.gene_name = gene
            msa_nuc = Readmsa.apply_dat_info_on_msa(msa_nuc, self.dat)
            self.logger.debug(f"{msa_nuc}")

        if "gen" in filetype and "nuc" in filetype:
            # remove some gen not included in nuc
            diff_name = list(set(msa_gen.get_sequence_names()) - set(msa_nuc.get_sequence_names()))
            if diff_name:
                self.logger.warning(f"Remove alleles doesn't exist in gen and nuc either: {diff_name}")
            msa_gen = msa_gen.remove(diff_name)

            # merge
            msa_merged = msa_gen.merge_exon(msa_nuc)
            self.logger.debug(f"{msa_merged}")
            return msa_merged
        elif "gen" in filetype:
            return msa_gen
        elif "nuc" in filetype:
            return msa_nuc
        else:
            return None


class KIRmsa(Familymsa):
    """
    A KIR interface that read MSA from MSF and KIR.dat

    Attributes:
        genes (dict of str, Genemsa): The msa object for each gene
    """

    def __init__(self, genes=[], filetype=["gen", "nuc"],
                 ipd_folder="KIR_v2100", version="2100"):
        super().__init__(genes, filetype, db_folder=ipd_folder, version=version)

    def _download_db(self, version="2100"):
        """
        Download the KIR to `IPDKIR`
        """
        # TODO: Auto find the latest version
        self._run_shell("git", "clone", "https://github.com/ANHIG/IPDKIR", self.db_folder)
        self._run_shell("git", "checkout", version, cwd=self.db_folder)

    def _get_name(self, search_name: str):
        """ Extract name from file pattern """
        arr_files = glob(search_name)
        return set([f.split("/")[-1].split("_")[0] for f in arr_files])

    def _list_db_gene(self, filetype) -> List[str]:
        """ List the gene in folder """
        if "gen" in filetype:
            names = names_gen = self._get_name(f"{self.db_folder}/msf/*_gen.msf")
        if "nuc" in filetype:
            # Most's of E nuc is different from dat record
            names = names_nuc = self._get_name(f"{self.db_folder}/msf/*_nuc.msf")
        if "gen" in filetype and "nuc" in filetype:
            names = names_gen & names_nuc
        return sorted(names)

    def _read_db_gene(self, gene: str, filetype: List[str]):
        """
        Read `{gene}_{filetype}.txt`.

        If both `gen` and `nuc` are given, it will merge them.

        The alignment data is stored in `self.genes[gene]` in `Genemsa` instance
        """
        if not hasattr(self, "dat"):
            self.logger.debug(f"Reading kir.dat")
            if os.path.exists(f"{self.db_folder}/KIR.dat"):
                self.dat = Readmsa.read_dat_block(f"{self.db_folder}/KIR.dat")
            else:
                self.dat = Readmsa.read_dat_block(f"{self.db_folder}/kir.dat")

        if "gen" in filetype:
            msa_gen = Readmsa.from_MSF_file(f"{self.db_folder}/msf/{gene}_gen.msf")
            msa_gen.seq_type = "gen"
            msa_gen.gene_name = gene
            msa_gen = Readmsa.apply_dat_info_on_msa(msa_gen, self.dat)
            self.logger.debug(f"Gen {msa_gen}")

        if "nuc" in filetype:
            msa_nuc = Genemsa(gene, seq_type="nuc")
            msa_nuc = Readmsa.from_MSF_file(f"{self.db_folder}/msf/{gene}_nuc.msf")
            msa_nuc.seq_type = "nuc"
            msa_nuc.gene_name = gene
            msa_nuc = Readmsa.apply_dat_info_on_msa(msa_nuc, self.dat)
            self.logger.debug(f"Nuc {msa_nuc}")

        if "gen" in filetype and "nuc" in filetype:
            # remove some gen not included in nuc
            diff_name = list(set(msa_gen.get_sequence_names()) - set(msa_nuc.get_sequence_names()))
            if diff_name:
                self.logger.warning(f"Remove alleles doesn't exist in gen and nuc either: {diff_name}")
            msa_gen = msa_gen.remove(diff_name)

            # specical case
            # exon 3 is pseudo exon
            # so, fill with gene's exon3
            gene_has_pseudo_exon3 = ["KIR2DL1", "KIR2DL2", "KIR2DL3", "KIR2DP1",
                                     "KIR2DS1", "KIR2DS2", "KIR2DS3", "KIR2DS4", "KIR2DS5"]
            if gene in gene_has_pseudo_exon3:

                exon3 = msa_gen.select_block([5])
                for name in set(msa_nuc.get_sequence_names()) - set(msa_gen.get_sequence_names()):
                    exon3.append(name, "-" * exon3.get_length())
                msas = msa_nuc.split()
                msa_nuc = msa_nuc.select_block(list(range(0, 2))) \
                          + exon3 \
                          + msa_nuc.select_block(list(range(3, len(msa_nuc.blocks))))
                msa_nuc.seq_type = "nuc"

            # merge
            msa_merged = msa_gen.merge_exon(msa_nuc)
            self.logger.debug(f"Merged {msa_merged}")
            return msa_merged

        elif "gen" in filetype:
            return msa_gen
        elif "nuc" in filetype:
            return msa_nuc
        else:
            return None
