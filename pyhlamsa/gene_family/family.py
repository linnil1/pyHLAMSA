import os
import logging
import subprocess
from typing import Any, List, Union, Set, Iterable

from ..gene import Genemsa, BlockInfo
from .. import msaio


GeneSet = Union[str, Iterable[str], None]
TypeSet = Union[str, Iterable[str]]


class Familymsa:
    """
    A abstract class to handle Gene Family

    When init, it will

    * Init parameter
    * Download db
    * Read db to Genemsa

    Attributes:
      genes (dict[str, Genemsa]):
        The dictionary use gene_name as key and msa object as value
    """

    def __init__(self, genes: GeneSet = None,
                 filetype: TypeSet = ["gen", "nuc"],
                 db_folder="dbpath", version="latest"):
        self.logger = logging.getLogger(__name__)
        self.db_folder = db_folder
        self.genes = {}

        # Download if not exist
        self._download(version)

        if genes is None:
            genes = self.list_db_gene(filetype)
        if isinstance(genes, str):
            genes = [genes]
        if isinstance(filetype, str):
            filetype = [filetype]
        filetype = set(filetype)

        # main
        for gene_name in genes:
            self.logger.info(f"Reading {gene_name}'s sequences")
            self.genes[gene_name] = self.read_db_gene(gene_name, filetype)

    def list_genes(self) -> List[str]:
        """ List all the gene's name in this family """
        return list(self.genes.keys())

    def __getitem__(self, index: str) -> Genemsa:
        """ Get specific gene's msa """
        return self.genes[index]

    def __iter__(self):
        """ Iter gene name like iter(dict) """
        return iter(self.genes)

    def items(self):
        """ list gene name and msa like dict.items() """
        return self.genes.items()

    def _download(self, version: str):
        """ Check before running download """
        if not os.path.exists(self.db_folder):
            self._download_db(version=version)
            if not os.path.exists(self.db_folder):
                raise ValueError("Fail to download")
        else:
            self.logger.info(f"{self.db_folder} exists")

    def _download_db(self, version: str):
        """ Abstract method: code for downloading your db """
        raise NotImplementedError

    def _run_shell(self, *args, cwd=None):
        """ Run shell code """
        self.logger.debug("Run " + " ".join(args))
        subprocess.run(args, cwd=cwd, check=True)

    def list_db_gene(self, filetype: TypeSet) -> List[str]:
        """ Abstract method: code for listing gene names """
        raise NotImplementedError

    def read_db_gene(self, gene: str, filetype: TypeSet) -> Genemsa:
        """ Abstract method: code for reading and merging function """
        raise NotImplementedError
