import os
from glob import glob
from .family import GeneSet, TypeSet, Familymsa
from ..gene import Genemsa
from ..utils import dat


class KIRmsa(Familymsa):
    """
    A KIR interface that read MSA from MSF and KIR.dat

    Attributes:
        genes (dict[str, Genemsa]): The msa object for each gene
    """

    def __init__(
        self,
        genes: GeneSet = None,
        filetype: TypeSet = ["gen", "nuc"],
        ipd_folder: str = "",
        version: str = "Latest",
    ):
        """
        Args:
            genes (str | list[str]): A list of genes you want to read.

                Set None if you want read all gene in HLA

            filetype (str | list[str] | set[str]): A list of filetype.

                If both `gen` and `nuc` are given, it will merge them automatically.

            ipd_folder (str): Path to your IPD/KIR folder

                You can manually download <https://github.com/ANHIG/IPDKIR> and
                checkout to specific branch

                Or it will automatically download the database with assigned `version` to
                `ipd_folder`. Default is `./kIR_v{verion}`

            version (str): IMGT version you want to download (2110 for version 2.11.0)

                If `ipd_folder` is existed, this value will be ignored.

                version='Latest' to get latest version, however sometime it cannot work
                because database may change the format, or contains bugs
                (e.g. 2.11.0 https://github.com/ANHIG/IPDKIR/issues/44)
        """
        # Why not version 2110 -> 2DL4,2DL5 has exon 4
        if not ipd_folder:
            ipd_folder = f"KIR_v{version}"
        super().__init__(genes, filetype, db_folder=ipd_folder, version=version)

    def _download_db(self, version: str = "Latest") -> None:
        """
        Download the KIR to `IPDKIR`
        """
        self._run_shell(
            "git",
            "clone",
            "--branch",
            version,
            "--single-branch",
            "https://github.com/ANHIG/IPDKIR",
            self.db_folder,
        )

    def _get_name(self, search_name: str) -> set[str]:
        """Extract name from file pattern"""
        arr_files = glob(search_name)
        return set([f.split("/")[-1].split("_")[0] for f in arr_files])

    def list_db_gene(self, filetype: TypeSet) -> list[str]:
        """List the gene in folder"""
        if "gen" in filetype:
            names = names_gen = self._get_name(f"{self.db_folder}/msf/*_gen.msf")
        if "nuc" in filetype:
            # Most's of E nuc is different from dat record
            names = names_nuc = self._get_name(f"{self.db_folder}/msf/*_nuc.msf")
        if "gen" in filetype and "nuc" in filetype:
            names = names_gen & names_nuc
        return sorted(names)

    def read_db_gene(self, gene: str, filetype: TypeSet) -> Genemsa:
        """
        Read `{gene}_{filetype}.msf and kir.dat`.

        If both `gen` and `nuc` are given, it will merge them.
        """
        if not hasattr(self, "dat"):
            self.logger.debug(f"Reading kir.dat")
            if os.path.exists(f"{self.db_folder}/KIR.dat"):
                self.dat = dat.read_dat_block(f"{self.db_folder}/KIR.dat")
            else:
                self.dat = dat.read_dat_block(f"{self.db_folder}/kir.dat")

        if "gen" in filetype:
            msa_gen = Genemsa.read_msf_file(f"{self.db_folder}/msf/{gene}_gen.msf")
            msa_gen.gene_name = gene
            msa_gen = dat.apply_dat_info_on_msa(msa_gen, self.dat, seq_type="gen")
            self.logger.debug(f"Gen {msa_gen}")

        if "nuc" in filetype:
            msa_nuc = Genemsa.read_msf_file(f"{self.db_folder}/msf/{gene}_nuc.msf")
            msa_nuc.gene_name = gene
            msa_nuc = dat.apply_dat_info_on_msa(msa_nuc, self.dat, seq_type="nuc")
            self.logger.debug(f"Nuc {msa_nuc}")

        if "gen" in filetype and "nuc" in filetype:
            # remove some gen not included in nuc
            diff_name = list(
                set(msa_gen.get_sequence_names()) - set(msa_nuc.get_sequence_names())
            )
            if diff_name:
                self.logger.warning(
                    f"Remove alleles doesn't exist in gen and nuc either: {diff_name}"
                )
            msa_gen = msa_gen.remove_allele(diff_name)

            # specical case
            # exon 3 is pseudo exon
            # so, fill with gene's exon3
            gene_has_pseudo_exon3 = [
                "KIR2DL1",
                "KIR2DL2",
                "KIR2DL3",
                "KIR2DP1",
                "KIR2DS1",
                "KIR2DS2",
                "KIR2DS3",
                "KIR2DS4",
                "KIR2DS5",
            ]
            if gene in gene_has_pseudo_exon3:

                exon3 = msa_gen.select_block([5])
                for name in set(msa_nuc.get_sequence_names()) - set(
                    msa_gen.get_sequence_names()
                ):
                    exon3.append(name, "-" * exon3.get_length())
                msa_nuc = (
                    msa_nuc.select_block(list(range(0, 2)))
                    + exon3
                    + msa_nuc.select_block(list(range(3, len(msa_nuc.blocks))))
                )

            # merge
            msa_merged = msa_gen.merge_exon(msa_nuc)
            self.logger.debug(f"Merged {msa_merged}")
            return msa_merged

        elif "gen" in filetype:
            return msa_gen
        elif "nuc" in filetype:
            return msa_nuc
        raise ValueError("gen or nuc are not exist in filetype")
