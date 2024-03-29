from glob import glob
from tempfile import TemporaryDirectory

from .family import GeneSet, TypeSet, Familymsa
from ..gene import Genemsa


class HLAmsa(Familymsa):
    """
    A HLA interface.

    This module read the HLA MSA from alignments/*.txt files,
    which use `|` to separate the intron/exon/UTR regions, but
    the labels are assumed with the order:
    `5UTR exon1 intron1 exon2 ... exonN 3UTR`

    Attributes:
      genes (dict[str, Genemsa]):
        The dictionary use gene_name as key and msa object as value
    """

    def __init__(
        self,
        genes: GeneSet = None,
        filetype: TypeSet = ["gen", "nuc"],
        imgt_alignment_folder: str = "",
        version: str = "Latest",
    ):
        """
        Args:
            genes (str | list[str]): A list of genes you want to read.

                Set None if you want read all gene in HLA

            filetype (str | list[str] | set[str]): A list of filetype.

                If both `gen` and `nuc` are given, it will merge them automatically.

            imgt_alignment_folder (str): Path to your IMGT-formatted MSA alignment file

                You can manually download the folder from
                <https://github.com/ANHIG/IMGTHLA/tree/Latest/alignments>
                or <http://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/>
                and unzip Alignments_Rel_3470.zip

                Otherwise it will automatically download the database to
                `imgt_alignment_folder`. Default is `./alignment_v{verion}`

            version (str): IMGT version you want to download (e.g. 3470 for 3.47.0)

                If `imgt_alignment_folder` is existed, this value will be ignored.
                Use `Latest` to get latest version.
        """
        if not imgt_alignment_folder:
            imgt_alignment_folder = f"alignment_v{version}"
        super().__init__(
            genes, filetype, db_folder=imgt_alignment_folder, version=version
        )

    def _download_db(self, version: str = "Latest") -> None:
        """
        Download the IMGTHLA alignments folder to `db_folder`

        I get the alignments folder from git instead of
        <http://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/>
        for better version controlling
        """
        with TemporaryDirectory() as tmp_dir:
            self._run_shell(
                "git",
                "clone",
                "--branch",
                version,
                "--single-branch",
                "https://github.com/ANHIG/IMGTHLA.git",
                tmp_dir,
            )
            self._run_shell("mv", tmp_dir + "/alignments", self.db_folder)

    def _get_name(self, search_name: str) -> set[str]:
        """Handy function to list names from file pattern"""
        arr_files = glob(search_name)
        return set([f.split("/")[-1].split("_")[0] for f in arr_files])

    def list_db_gene(self, filetype: TypeSet) -> list[str]:
        """List the gene in folder"""
        drb = set(["DRB1", "DRB3", "DRB4", "DRB5"])
        if "gen" in filetype:
            names = names_gen = self._get_name(f"{self.db_folder}/*_gen.txt") | drb
        if "nuc" in filetype:
            names = names_nuc = self._get_name(f"{self.db_folder}/*_nuc.txt") | drb
        if "gen" in filetype and "nuc" in filetype:
            names = names_gen & names_nuc
            # * E exon7 is ambiguous, exon8 is gone (Also exon8 is not pseudo exon)
            #    <E gen alleles=258 block=5UTR(301) exon1(64) intron1(130)
            #        exon2(270) intron2(244) exon3(276) intron3(621) exon4(276)
            #        intron4(124) exon5(117) intron5(751) exon6(33) intron6(104)
            #        exon7(43) intron7(165) exon8(5) 3UTR(356)>
            #    <E nuc alleles=262 block=exon1(64) exon2(270) exon3(276)
            #        exon4(276) exon5(117) exon6(33) exon7(41)>
            if "gen" in filetype and "nuc" in filetype:
                names = names - set(["E"])
        return list(sorted(names))

    def read_db_gene(self, gene: str, filetype: TypeSet) -> Genemsa:
        """
        Read `{gene}_{filetype}.txt`.

        If both `gen` and `nuc` are given, it will merge them.
        """
        if "gen" in filetype:
            msa_gen = Genemsa.read_imgt_alignment(f"{self.db_folder}/{gene}_gen.txt")
            msa_gen.gene_name = gene
            if gene != "P":
                # P: special case: has even block
                # <P gen alleles=5 block=(475) (261) (589) (276) (124) (117)
                #                        (412) (33) (150) (48) (163) (297)>
                msa_gen.assume_label("gen")
            self.logger.debug(f"Gen {msa_gen}")
        if "nuc" in filetype:
            # Special Case: DRB* nuc are in DRB_nuc.txt
            if gene.startswith("DRB"):
                msa_nuc = Genemsa.read_imgt_alignment(f"{self.db_folder}/DRB_nuc.txt")
                msa_nuc = msa_nuc.select_allele(gene + ".*")
            else:
                msa_nuc = Genemsa.read_imgt_alignment(
                    f"{self.db_folder}/{gene}_nuc.txt"
                )

            msa_nuc.gene_name = gene
            msa_nuc.assume_label("nuc")
            self.logger.debug(f"Nuc {msa_nuc}")

        if "gen" in filetype and "nuc" in filetype:
            # remove some gene
            diff_name = list(
                set(msa_gen.get_sequence_names()) - set(msa_nuc.get_sequence_names())
            )
            if diff_name:
                self.logger.warning(
                    f"Remove alleles doesn't exist in gen and nuc either: {diff_name}"
                )
            msa_gen = msa_gen.remove_allele(diff_name)

            # merge
            msa_merged = msa_gen.merge_exon(msa_nuc)
            self.logger.debug(f"Merged {msa_merged}")
            return msa_merged

        elif "gen" in filetype:
            return msa_gen
        elif "nuc" in filetype:
            return msa_nuc
        raise ValueError("gen or nuc are not exist in filetype")
