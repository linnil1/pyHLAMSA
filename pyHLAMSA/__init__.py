from __future__ import annotations
import os
import logging
import subprocess
from glob import glob
from typing import List, Set
from Bio import SeqIO

from .Genemsa import Genemsa, BlockInfo, IndexInfo
from . import Readmsa


class Familymsa:
    """
    A abstract class to handle Gene Family

    When init, it will

    * Init parameter
    * Download db
    * Read db to Genemsa

    Attributes:
        genes (dict): The dictionary use gene_name as key and msa object as value
    """

    def __init__(self, genes=[], filetype=["gen", "nuc"],
                 db_folder="alignments", version="latest"):
        self.logger = logging.getLogger(__name__)
        self.db_folder = db_folder
        self.genes = {}

        # Download if not exist
        self._download(version)

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
        for gene_name in genes:
            self.logger.info(f"Reading {gene_name}'s sequences")
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
        genes (dict): The dictionary use gene_name as key and msa object as value
    """
    def __init__(self, genes=[], filetype=["gen", "nuc"],
                 imgt_alignment_folder=None, version="3470"):
        """
        Args:
            genes (str or list of str): A list of genes you want to read.

                Leave Empty if you want read all gene in HLA

            filetype (str or list of str): A list of filetype.

                If both `gen` and `nuc` are given, it will merge them automatically.

            imgt_alignment_folder (str): Path to your IMGT-formatted MSA alignment file

                You can manually download from
                <http://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/>
                and unzip Alignments_Rel_3470.zip

                Or it will automatically download the database with assigned `version` to
                `imgt_alignment_folder`. Default is `./alignment_v{verion}`

            version (str): IMGT version you want to download

                If `imgt_alignment_folder` is existed, this value will be ignored.
        """
        if imgt_alignment_folder is None:
            imgt_alignment_folder = f"alignment_v{version}"
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
            # * E exon7 is ambiguous, exon8 is gone (Also exon8 is not pseudo exon)
            #    <E gen alleles=258 block=5UTR(301) exon1(64) intron1(130)
            #        exon2(270) intron2(244) exon3(276) intron3(621) exon4(276)
            #        intron4(124) exon5(117) intron5(751) exon6(33) intron6(104)
            #        exon7(43) intron7(165) exon8(5) 3UTR(356)>
            #    <E nuc alleles=262 block=exon1(64) exon2(270) exon3(276)
            #        exon4(276) exon5(117) exon6(33) exon7(41)>
            if "gen" in filetype and "nuc" in filetype:
                names = names - set(["E"])
        return list(sorted(names))[10:]

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
                # <P gen alleles=5 block=(475) (261) (589) (276) (124) (117)
                #                        (412) (33) (150) (48) (163) (297)>
                msa_gen._assume_label()
            self.logger.debug(f"Gen {msa_gen}")
        if "nuc" in filetype:
            # Special Case: DRB* nuc are in DRB_nuc.txt
            if gene.startswith("DRB"):
                msa_nuc = Readmsa.from_alignment_file(f"{self.db_folder}/DRB_nuc.txt")
                msa_nuc = msa_nuc.select_allele(gene + ".*")
            else:
                msa_nuc = Readmsa.from_alignment_file(f"{self.db_folder}/{gene}_nuc.txt")

            msa_nuc.seq_type = "nuc"
            msa_nuc.gene_name = gene
            msa_nuc._assume_label()
            self.logger.debug(f"Nuc {msa_nuc}")

        if "gen" in filetype and "nuc" in filetype:
            # remove some gene
            diff_name = list(set(msa_gen.get_sequence_names())
                             - set(msa_nuc.get_sequence_names()))
            if diff_name:
                self.logger.warning(
                    f"Remove alleles doesn't exist in gen and nuc either: {diff_name}")
            msa_gen = msa_gen.remove(diff_name)

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


class HLAmsaEX(Familymsa):
    """
    A HLA interface but read gene from MSF and
    read intron/exon information from hla.dat

    I think this one is more reliable

    Attributes:
        genes (dict): The dictionary use gene_name as key and msa object as value
    """

    def __init__(self, genes=[], filetype=["gen", "nuc"],
                 imgt_folder=None, version="3470"):
        """
        Args:
            genes (str or list of str): A list of genes you want to read.

                Leave Empty if you want read all gene in HLA

            filetype (str or list of str): A list of filetype.

                If both `gen` and `nuc` are given, it will merge them automatically.

            imgt_folder (str): Path to your IMGT/HLA root folder

                You can manually download <https://github.com/ANHIG/IMGTHLA/> and
                checkout to specific branch

                Or it will automatically download the database with assigned `version` to
                `imgt_folder`. Default is `./IMGT_v{verion}`

            version (str): IMGT version you want to download

                If `imgt_folder` is existed, this value will be ignored.

                You can use `Latest` to get the latest version
        """
        if imgt_folder is None:
            imgt_folder = f"IMGT_v{version}"
        super().__init__(genes, filetype,
                         db_folder=imgt_folder, version=version)

    def _download_db(self, download=True, version="3430"):
        """
        Download the IMGT/HLA msf and hla.dat to folder `imgt_folder`
        """
        # TODO: Auto find the latest version
        self._run_shell("git", "clone",
                        "https://github.com/ANHIG/IMGTHLA.git", self.db_folder)
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
            # HLA-E nuc doessn't have exon8
            names = names_nuc = (self._get_name(f"{self.db_folder}/msf/*_nuc.msf") | DRB
                                 - set(["E"]))
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
            diff_name = list(set(msa_gen.get_sequence_names())
                             - set(msa_nuc.get_sequence_names()))
            if diff_name:
                self.logger.warning(
                    f"Remove alleles doesn't exist in gen and nuc either: {diff_name}")
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
                 ipd_folder=None, version="2100"):
        """
        Args:
            genes (str or list of str): A list of genes you want to read.

                Leave Empty if you want read all gene in HLA

            filetype (str or list of str): A list of filetype.

                If both `gen` and `nuc` are given, it will merge them automatically.

            ipd_folder (str): Path to your IPD/KIR folder

                You can manually download <https://github.com/ANHIG/IPDKIR> and
                checkout to specific branch

                Or it will automatically download the database with assigned `version` to
                `ipd_folder`. Default is `./kIR_v{verion}`

            version (str): IMGT version you want to download

                If `ipd_folder` is existed, this value will be ignored.

                version='Latest' to get latest version, however sometime it cannot work
                because database may change the format, or contains bugs
        """
        # Why not version 2110 -> 2DL4,2DL5 has exon 4
        if ipd_folder is None:
            ipd_folder = f"KIR_v{version}"
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
            msa_nuc = Readmsa.from_MSF_file(f"{self.db_folder}/msf/{gene}_nuc.msf")
            msa_nuc.seq_type = "nuc"
            msa_nuc.gene_name = gene
            msa_nuc = Readmsa.apply_dat_info_on_msa(msa_nuc, self.dat)
            self.logger.debug(f"Nuc {msa_nuc}")

        if "gen" in filetype and "nuc" in filetype:
            # remove some gen not included in nuc
            diff_name = list(set(msa_gen.get_sequence_names()) -
                             set(msa_nuc.get_sequence_names()))
            if diff_name:
                self.logger.warning(
                    f"Remove alleles doesn't exist in gen and nuc either: {diff_name}")
            msa_gen = msa_gen.remove(diff_name)

            # specical case
            # exon 3 is pseudo exon
            # so, fill with gene's exon3
            gene_has_pseudo_exon3 = ["KIR2DL1", "KIR2DL2", "KIR2DL3", "KIR2DP1",
                                     "KIR2DS1", "KIR2DS2", "KIR2DS3", "KIR2DS4",
                                     "KIR2DS5"]
            if gene in gene_has_pseudo_exon3:

                exon3 = msa_gen.select_block([5])
                for name in (set(msa_nuc.get_sequence_names())
                             - set(msa_gen.get_sequence_names())):
                    exon3.append(name, "-" * exon3.get_length())
                msas = msa_nuc.split()
                msa_nuc = (msa_nuc.select_block(list(range(0, 2)))
                           + exon3
                           + msa_nuc.select_block(list(range(3, len(msa_nuc.blocks)))))
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


class CYPmsa(Familymsa):
    """
    A CYP interface that read MSA from fasta and haplotypes.tsv

    Attributes:
        genes (dict of str, Genemsa): The msa object for each gene
    """

    def __init__(self, genes=[], filetype=[],
                 pharmvar_folder=None, version="5.1.10"):
        """
        Args:
            genes (str or list of str): A list of genes you want to read.

                Leave Empty if you want read all gene in HLA

            filetype: Ignore

            pharmvar_folder (str): Path to your pharmvar folder

                You should manually download <https://www.pharmvar.org/download>
                and unzip it.

            version (str): PharmVar version you want to download

                Not works now
        """
        if pharmvar_folder is None:
            pharmvar_folder = f"pharmvar-{version}"
        super().__init__(genes, filetype, db_folder=pharmvar_folder, version=version)

    def _download_db(self, version):
        """
        Download the CYP from https://www.pharmvar.org/download
        """
        raise ValueError("You should download CYP genes from "
                         "https://www.pharmvar.org/download")

    def _list_db_gene(self, filetype) -> List[str]:
        """ List the gene in folder """
        return sorted(os.listdir(self.db_folder))

    def _read_db_gene(self, gene: str, filetype: List[str]):
        """
        Read `.haplotypes.fasta` and `.haplotypes.tsv`

        `filetype` will be ignored now
        """
        # read fasta
        ref_seqs = {}
        for seq in SeqIO.parse(f"{self.db_folder}/{gene}/{gene}.haplotypes.fasta",
                               "fasta"):
            # split allele name
            # e.g.  rs75017182, rs56038477 PV01077 NG_008807.2 PharmVar Version:5.1.10
            for name in seq.description.replace(', ', ',').split(' ')[0].split(','):
                ref_seqs[name.strip()] = str(seq.seq)
        self.logger.debug(f"Read sequence {ref_seqs.keys()}")

        # read tsv
        # split by '\t' and ignore header
        table = open(glob(f"{self.db_folder}/{gene}/RefSeqGene/*.haplotypes.tsv")[0])
        table = [i.strip().split('\t') for i in table if not i.startswith("#")][1:]

        # Get Reference sequence
        reference_seq = None
        reference_name = None
        for i in table:
            if i[3] == "REFERENCE" and ref_seqs.get(i[0])[0]:
                if reference_seq:
                    assert reference_seq == ref_seqs[i[0]]
                    continue
                reference_name = i[0]
                reference_seq = ref_seqs[reference_name]
        if not reference_seq:
            # hard-coded for DPYD, because it's reference `NG_008807.2` is missing
            if gene == "DPYD":
                reference_name = "rs111858276"
                reference_seq = ref_seqs[reference_name]
                reference_seq = reference_seq[:376460 - 1] + "A" + reference_seq[376460:]
            else:
                raise ValueError(f"Not reference found in {gene}")

        # Fill with Reference and reference is the first one
        alleles = {reference_name: None}
        alleles.update({
            allele_name: reference_seq
            for allele_name in set(i[0] for i in table)
        })
        self.logger.debug(f"Read vcf {alleles.keys()}")

        # Remove duplciate
        # in 5.1.10, there exist two same row
        # CYP2D6*149  CYP2D6  rs59421388  NG_008376.4 8203    8203    G   A   substitution
        # CYP2D6*149  CYP2D6  rs59421388  NG_008376.4 8203    8203    G   A   substitution
        table = list(set(map(tuple, filter(lambda i: i[3] != "REFERENCE", table))))

        # Reconstruct msa from VCF variant
        # VCF type: substitution and deleteion, insertion
        # Table Header: ['Haplotype Name', 'Gene', 'rsID', 'ReferenceSequence',
        #   'Variant Start', 'Variant Stop', 'Reference Allele', 'Variant Allele', 'Type']
        # Insertion will move the index, so insert from last position first
        for i in sorted(map(list, table), key=lambda i: -int(i[4])):
            if i[8] == "substitution":
                pos = int(i[4]) - 1
                end = pos + 1
            elif i[8] == "deletion":
                pos = int(i[4]) - 1
                end = pos + len(i[6])
                i[7] = '-' * len(i[6])
            elif i[8] == "insertion":
                pos = int(i[4])
                end = int(i[4]) + len(i[7])
                i[6] = '-' * len(i[7])
                # check the gap is already inserted by others
                # if false: insert '-' for all alleles
                # else: replace the '-' with insertion seq
                if set(alleles[i[0]][pos:end]) != set("-"):
                    for allele_name in alleles:
                        alleles[allele_name] = alleles[allele_name][:pos] + i[6] \
                                               + alleles[allele_name][pos:]
            else:
                raise ValueError(f"Unknown Type {i[8]}")

            assert alleles[i[0]][pos:end] == i[6]
            alleles[i[0]] = alleles[i[0]][:pos] + i[7] + alleles[i[0]][end:]

        # split allele name
        # e.g. rs75017182, rs56038477  DPYD    rs75017182  NG_008807.2 346167
        for allele_name in list(alleles.keys()):
            if ', ' in allele_name:
                for an in allele_name.split(","):
                    alleles[an.strip()] = alleles[allele_name]
                del alleles[allele_name]

        # double check
        length = None
        for allele_name in alleles:
            if length is None:
                length = len(alleles[allele_name])
            else:
                assert length == len(alleles[allele_name])
            if alleles[allele_name].replace("-", "") != ref_seqs[allele_name]:
                raise ValueError(f"{allele_name} is not same as reference")

        msa = Genemsa(gene, "gen")
        msa.alleles = alleles
        msa.blocks = [BlockInfo(length=length, type="gene", name="gene")]
        return msa.reset_index()
