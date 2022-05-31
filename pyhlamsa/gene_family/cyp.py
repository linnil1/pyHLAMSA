import os
from glob import glob
from typing import List, cast
from Bio import SeqIO

from .family import Familymsa, TypeSet, GeneSet, Genemsa, BlockInfo


class CYPmsa(Familymsa):
    """
    A CYP interface that read MSA from fasta and haplotypes.tsv

    Attributes:
        genes (dict[str, Genemsa]): The msa object for each gene
    """

    def __init__(self, genes: GeneSet = None,
                 filetype: TypeSet = [],
                 pharmvar_folder="", version="5.1.10"):
        """
        Args:
            genes (str | list[str]): A list of genes you want to read.

                Set None if you want read all gene in HLA

            filetype: Ignore

            pharmvar_folder (str): Path to your pharmvar folder

                You should manually download <https://www.pharmvar.org/download>
                and unzip it.

            version (str): PharmVar version you want to download

                Not works now
        """
        if not pharmvar_folder:
            pharmvar_folder = f"pharmvar-{version}"
        super().__init__(genes, filetype, db_folder=pharmvar_folder, version=version)

    def _download_db(self, version=""):
        """
        Download the CYP from https://www.pharmvar.org/download
        """
        raise ValueError("You should download CYP genes from "
                         "https://www.pharmvar.org/download")

    def list_db_gene(self, filetype: TypeSet = []) -> List[str]:
        """ List the gene in folder """
        return sorted(os.listdir(self.db_folder))

    def read_db_gene(self, gene: str, filetype: TypeSet = []) -> Genemsa:
        """
        Read `{gene}/{gene}.haplotypes.fasta` and `{gene}/{gene}.haplotypes.tsv`

        `filetype` will be ignored now
        """
        # read fasta
        ref_seqs = {}
        for seq in SeqIO.parse(f"{self.db_folder}/{gene}/{gene}.haplotypes.fasta",
                               "fasta"):
            # split allele name
            # e.g.  rs75017182, rs56038477 PV01077 NG_008807.2 PharmVar Version:5.1.10
            for name in seq.description.replace(", ", ",").split(" ")[0].split(","):
                ref_seqs[name.strip()] = str(seq.seq)
        self.logger.debug(f"Read sequence {ref_seqs.keys()}")

        # read tsv
        # split by '\t' and ignore header
        var_text = open(glob(f"{self.db_folder}/{gene}/RefSeqGene/*.haplotypes.tsv")[0])
        table = [i.strip().split("\t") for i in var_text if not i.startswith("#")][1:]

        # Get Reference sequence
        reference_seq = None
        reference_name = None
        for i in table:
            if i[0] not in ref_seqs:
                continue
            if i[3] == "REFERENCE" and ref_seqs[i[0]][0]:
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

        assert reference_name is not None
        # Fill with Reference and reference is the first one
        alleles = {reference_name: ""}
        alleles.update({
            allele_name: reference_seq
            for allele_name in set(i[0] for i in table)
        })
        self.logger.debug(f"Read vcf {alleles.keys()}")

        # Remove duplciate
        # in 5.1.10, there exist two same row
        # CYP2D6*149  CYP2D6  rs59421388  NG_008376.4 8203    8203    G   A   substitution
        # CYP2D6*149  CYP2D6  rs59421388  NG_008376.4 8203    8203    G   A   substitution
        table_unique1 = filter(lambda i: i[3] != "REFERENCE", table)
        table_unique2 = cast(tuple, map(tuple, table_unique1))  # list to tuple
        table_unique = cast(list, map(list, set(list(table_unique2))))  # tuple to list

        # Reconstruct msa from VCF variant
        # VCF type: substitution and deleteion, insertion
        # Table Header: ['Haplotype Name', 'Gene', 'rsID', 'ReferenceSequence',
        #   'Variant Start', 'Variant Stop', 'Reference Allele', 'Variant Allele', 'Type']
        # Insertion will move the index, so insert from last position first
        for i in sorted(table_unique, key=lambda i: -int(i[4])):
            if i[8] == "substitution":
                pos = int(i[4]) - 1
                end = pos + 1
            elif i[8] == "deletion":
                pos = int(i[4]) - 1
                end = pos + len(i[6])
                i[7] = "-" * len(i[6])
            elif i[8] == "insertion":
                pos = int(i[4])
                end = int(i[4]) + len(i[7])
                i[6] = "-" * len(i[7])
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

        # split allele name (alleles.key() will change)
        # e.g. rs75017182, rs56038477  DPYD    rs75017182  NG_008807.2 346167
        for allele_name in list(alleles.keys()):
            if ", " in allele_name:
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

        assert length
        msa = Genemsa(gene)
        msa.alleles = alleles
        msa.blocks = [BlockInfo(length=length)]
        return msa.assume_label("other")
