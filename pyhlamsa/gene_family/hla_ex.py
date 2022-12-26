from glob import glob
import subprocess

from .family import Familymsa, Genemsa, GeneSet, TypeSet
from ..utils import dat


class HLAmsaEX(Familymsa):
    """
    A HLA interface but read gene from MSF and
    read intron/exon information from hla.dat

    I think this one is more reliable

    Attributes:
      genes (dict[str, Genemsa]):
        The dictionary use gene_name as key and msa object as value
    """

    def __init__(
        self,
        genes: GeneSet = None,
        filetype: TypeSet = ["gen", "nuc"],
        imgt_folder: str = "",
        version: str = "Latest",
    ):
        """
        Args:
            genes (str | list[str]): A list of genes you want to read.

                Set None if you want read all gene in HLA

            filetype (str | list[str] | set[str]): A list of filetype.

                If both `gen` and `nuc` are given, it will merge them automatically.

            imgt_folder (str): Path to your IMGT/HLA root folder

                You can manually download <https://github.com/ANHIG/IMGTHLA/> and
                checkout to specific branch

                Or it will automatically download the database with assigned `version` to
                `imgt_folder`. Default is `./IMGT_v{verion}`

            version (str): IMGT version you want to download (e.g. 3470 for 3.47.0)

                If `imgt_folder` is existed, this value will be ignored.

                You can use `Latest` to get the latest version
        """
        if not imgt_folder:
            imgt_folder = f"IMGT_v{version}"
        super().__init__(genes, filetype, db_folder=imgt_folder, version=version)

    def _download_db(self, version: str = "Latest"):
        """
        Download the IMGT/HLA msf and hla.dat to folder `imgt_folder`
        """
        try:
            self._run_shell("git", "lfs")
        except subprocess.CalledProcessError:
            self.logger.error(f"git lfs not installed (Try apt install git-lfs)")
            exit()

        self._run_shell(
            "git",
            "clone",
            "--branch",
            version,
            "--single-branch",
            "https://github.com/ANHIG/IMGTHLA.git",
            self.db_folder,
        )
        self._run_shell("git", "lfs", "pull", cwd=self.db_folder)

    def _get_name(self, search_name: str) -> set[str]:
        """Extract name from file pattern"""
        arr_files = glob(search_name)
        return set([f.split("/")[-1].split("_")[0] for f in arr_files])

    def list_db_gene(self, filetype: TypeSet) -> list[str]:
        """List the gene in folder"""
        drb = set(["DRB1", "DRB3", "DRB4", "DRB5"])
        if "gen" in filetype:
            names = names_gen = self._get_name(f"{self.db_folder}/msf/*_gen.msf") | drb
        if "nuc" in filetype:
            names = names_nuc = self._get_name(f"{self.db_folder}/msf/*_nuc.msf") | drb
        if "gen" in filetype and "nuc" in filetype:
            # HLA-E nuc has 7 bases differences from exon7 and exon8
            names = (names_gen & names_nuc) - set("E")
        return sorted(names)

    def read_db_gene(self, gene: str, filetype: TypeSet) -> Genemsa:
        """
        Read `msf/{gene}_{filetype}.msf` and hla.dat.

        If both `gen` and `nuc` are given, it will merge them.
        """
        if not hasattr(self, "dat"):
            self.logger.debug(f"Reading hla.dat")
            self.dat = dat.read_dat_block(f"{self.db_folder}/hla.dat")

        if "gen" in filetype:
            msa_gen = Genemsa.read_msf_file(f"{self.db_folder}/msf/{gene}_gen.msf")
            msa_gen.gene_name = gene
            msa_gen = dat.apply_dat_info_on_msa(msa_gen, self.dat, seq_type="gen")
            self.logger.debug(f"{msa_gen}")

        if "nuc" in filetype:
            # Special Case: DRB* nuc are in DRB_nuc.txt
            if gene.startswith("DRB"):
                msa_nuc = Genemsa.read_msf_file(f"{self.db_folder}/msf/DRB_nuc.msf")
                msa_nuc = msa_nuc.select_allele(gene + ".*")
            else:
                msa_nuc = Genemsa.read_msf_file(f"{self.db_folder}/msf/{gene}_nuc.msf")
            msa_nuc.gene_name = gene
            msa_nuc = dat.apply_dat_info_on_msa(msa_nuc, self.dat, seq_type="nuc")
            self.logger.debug(f"{msa_nuc}")

        if "gen" in filetype and "nuc" in filetype:
            # remove some gen not included in nuc
            diff_name = list(
                set(msa_gen.get_sequence_names()) - set(msa_nuc.get_sequence_names())
            )
            if diff_name:
                self.logger.warning(
                    f"Remove alleles existed in gen but not in nuc: {diff_name}"
                )
            msa_gen = msa_gen.remove(diff_name)

            # merge
            msa_merged = msa_gen.merge_exon(msa_nuc)
            self.logger.debug(f"{msa_merged}")
            return msa_merged
        elif "gen" in filetype:
            return msa_gen
        elif "nuc" in filetype:
            return msa_nuc
        raise ValueError("gen or nuc are not exist in filetype")
