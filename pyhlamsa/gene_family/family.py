import os
import logging
import subprocess
from abc import ABC
from typing import Any, List

from ..gene import Genemsa, BlockInfo
from .. import msaio


class Familymsa(ABC):
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
        if isinstance(genes, str):
            genes = [genes]
        if isinstance(filetype, str):
            filetype = [filetype]

        filetype = set(filetype)
        if not genes:  # if empty -> read all
            genes = self._list_db_gene(filetype)

        # main
        for gene_name in genes:
            self.logger.info(f"Reading {gene_name}'s sequences")
            self.genes[gene_name] = self._read_db_gene(gene_name, filetype)

    def list_genes(self) -> List[str]:
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

    def _download_db(self, version: str):
        """ Abstract method: code for downloading your db """
        raise NotImplementedError

    def _run_shell(self, *args, cwd=None):
        """ Run shell code """
        self.logger.debug("Run " + " ".join(args))
        with subprocess.Popen(args, cwd=cwd) as proc:
            proc.wait()

    def _list_db_gene(self, filetype: Any) -> List[str]:
        """ Abstract method: code for listing gene names """
        raise NotImplementedError

    def _read_db_gene(self, gene, filetype: Any) -> Genemsa:
        """ Abstract method: code for reading and merging function """
        raise NotImplementedError
