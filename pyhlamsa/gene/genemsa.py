"""
Gene Module.

All the msa operation are merged into Genemsa
"""
from __future__ import annotations
from .base import GenemsaBase
from .block_operation import GenemsaBlockOp
from .column_operation import GenemsaColumnOp
from .type_convertion import GenemsaConverter
from .allele_operation import GenemsaAlleleOp
from .text_operation import GenemsaTextOp


class Genemsa(GenemsaAlleleOp,
              GenemsaConverter,
              GenemsaTextOp,
              GenemsaColumnOp, GenemsaBlockOp,
              GenemsaBase):
    """
    Genemsa: core object in pyhlamsa

    It can do many msa operation while
    automatically updating index (column position) and block(intron exon region).

    Basic usage:
    ``` python
    from pyhlamsa import Genemsa
    ```
    """
    def select_complete(self) -> Genemsa:
        """ Select non exon-only sequences (No `E` in the sequence)"""
        new_msa = self.copy(copy_allele=False)
        new_msa.alleles = {allele: seq for allele, seq in self.alleles.items()
                           if "E" not in seq}
        return new_msa

    def select_incomplete(self) -> Genemsa:
        """ Select exon-only sequences (`E` exists in the sequence)"""
        new_msa = self.copy(copy_allele=False)
        new_msa.alleles = {allele: seq for allele, seq in self.alleles.items()
                           if "E" in seq}
        return new_msa

    def fill_incomplete(self, ref_allele: str) -> Genemsa:
        """ Fill the `E` in exon-only sequences with ref_allele sequence (inplace) """
        if ref_allele not in self.alleles:
            raise KeyError(f"{ref_allele} not found")

        ref_seq = self.alleles[ref_allele]
        for allele, seq in self.alleles.items():
            if "E" in seq:
                self.alleles[allele] = "".join(
                    [seq[i] if seq[i] != "E" else ref_seq[i]
                     for i in range(len(seq))])
        return self

    def merge_exon(self, msa_nuc: Genemsa) -> Genemsa:
        """
        Merge nuc MSA into gen MSA

        It's allow that nuc MSA has new allele name than gen MSA,
        Genemsa will add the sequence in MSA, and the intron will fill by `E`

        If the exon part of gen MSA is differnet (e.g. less gapped) from nuc MSA,
        Genemsa will try to merge if it can

        Note that the index will be reset

        Example:
          ```
          # case1
          msa_gen:
            1: "AA|TT|CC",
            2: "AA|TC|CC",
          msa_nuc:
            3: "TC",
          After merge:
            1: "AA|TT|CC",
            2: "AA|TC|CC",
            3: "EE|TC|EE"
          ```

          ```
          # case2
          msa_gen:
            1: "AA|TT|CC",
            2: "AA|TC|CC",
          msa_nuc:
            1: "TT-",
            2: "T-C",
            4: "TTC",
          After merge:
            1: "AA|TT-|CC",
            2: "AA|T-C|CC",
            4: "EE|TTC|EE"
          ```
        """
        # A mapping from gen name to nuc index
        nuc_name_index = {b.name: i for i, b in enumerate(msa_nuc.blocks)
                          if b.type == "exon"}

        # check it's one-to-one mapping
        exon_set = set(b.name for b in self.blocks if b.type == "exon")
        if set(nuc_name_index.keys()) != exon_set:
            raise ValueError(f"Cannot match blocks: "
                             f"gen={exon_set} nuc={nuc_name_index.keys()}")

        # create new msa and make sure the order of alleles
        new_msa = Genemsa(self.gene_name, reference=self.reference)
        new_msa.alleles = {name: "" for name in self.get_sequence_names()}
        new_msa.alleles.update({name: "" for name in msa_nuc.get_sequence_names()})

        # allele names
        gen_names = set(self.get_sequence_names())
        nuc_names = set(msa_nuc.get_sequence_names())
        exclude_name = set()  # type: set[str]

        # block-wise
        msas_gen = self.split()
        msas_nuc = msa_nuc.split()
        for i_gen in range(len(self.blocks)):
            # intron -> fill with E
            if self.blocks[i_gen].name not in nuc_name_index:
                for name in nuc_names - gen_names:
                    msas_gen[i_gen].append(name,
                                           "E" * self.blocks[i_gen].length)
                new_msa += msas_gen[i_gen].remove(list(exclude_name))
            # exon -> check before merge
            else:
                i_nuc = nuc_name_index[self.blocks[i_gen].name]
                # content change or length change
                if (msas_nuc[i_nuc].get_length() != msas_gen[i_gen].get_length()
                    or any(msas_nuc[i_nuc].get(name) != msas_gen[i_gen].get(name)
                           for name in (nuc_names & gen_names))):
                    # check before merge
                    if len(gen_names - nuc_names):
                        raise ValueError(
                            f"Some alleles doesn't exist in nuc MSA: "
                            f"{gen_names - nuc_names}")

                    diff_name = filter(lambda name:
                                       msas_nuc[i_nuc].get(name).replace("-", "")
                                       != msas_gen[i_gen].get(name).replace("-", ""),
                                       gen_names)
                    diff_names = list(diff_name)
                    if diff_names:
                        self.logger.warning(
                            f"Some exon sequences in gen MSA "
                            f"is not same as in nuc MSA "
                            f"{self.blocks[i_gen].name}: {diff_names}")
                        new_msa.remove(diff_names)
                        exclude_name.update(diff_names)
                new_msa += msas_nuc[i_nuc].remove(list(exclude_name))
        return new_msa.reset_index()

    def assume_label(self, seq_type="gen") -> Genemsa:
        """
        It will automatically generate the block's label
        according on `seq_type`. (Inplace)

        seq_type:

          * gen: 5UTR-exon1-intron1-exon2-...-exon9-3UTR
          * nuc: exon1-exon2-...-exon9
          * other: block1-block2-block3-...
        """
        if seq_type == "gen":
            assert len(self.blocks) % 2 == 1 and len(self.blocks) >= 3
            for i in range(len(self.blocks)):
                if i % 2:
                    self.blocks[i].type = "exon"
                    self.blocks[i].name = f"exon{i // 2 + 1}"
                else:
                    self.blocks[i].type = "intron"
                    self.blocks[i].name = f"intron{i // 2}"

            self.blocks[0].type = "five_prime_UTR"
            self.blocks[0].name = "5UTR"
            self.blocks[-1].type = "three_prime_UTR"
            self.blocks[-1].name = "3UTR"
        elif seq_type == "nuc":
            for i in range(len(self.blocks)):
                self.blocks[i].type = "exon"
                self.blocks[i].name = f"exon{i+1}"
        else:
            for i in range(len(self.blocks)):
                self.blocks[i].type = "gene_fragment"
                self.blocks[i].name = f"block{i+1}"

        # inplace to reset the index
        self.index = self.reset_index().index
        return self
