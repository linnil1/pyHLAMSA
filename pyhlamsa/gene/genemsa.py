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

    It can do many msa operation while
    automatically updating index (column position) and block(intron exon region).

    Basic usage:
    ``` python
    from pyhlamsa import Genemsa
    ```
    """

    pass
