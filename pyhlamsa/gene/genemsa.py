"""
Gene Module.

All the msa operation are merged into Genemsa
"""
from __future__ import annotations
from .base import GenemsaBase
from .io_operation import GenemsaIO
from .text_operation import GenemsaTextOp


class Genemsa(
    GenemsaIO,
    GenemsaTextOp,
    GenemsaBase,
):
    """
    Genemsa: core object in pyhlamsa

    More details written in GenemsaBase

    Basic usage:
    ``` python
    from pyhlamsa import Genemsa
    msa = Genemsa("test1")
    msa.assume_label("gen", blocks_length=[3,5,1])
    msa.append("A", "GGAGGA--CA")
    msa.append("B", "GGGGGAA--A")
    msa.print_alignment()
                      1   4     9
                      |   |     |
    A                 GGA|GGA--|CA
    B                 GGG|GGAA-|-A
    ```
    """

    pass
