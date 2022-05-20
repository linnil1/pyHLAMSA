"""
Utility for cigar calculation
"""
from typing import List, Tuple


def calculate_cigar(ref: str, seq: str) -> List[Tuple[str, int]]:
    """
    Compare two sequences and output cigar

    Example:
      ```
      reference_seq:    AATTAT
      target_seq:       AC--AT
      cigar: [(M, 1), (X, 1), (D, 2), (M, 2)]
      ```
    """
    if len(ref) != len(seq):
        raise ValueError("Two sequences doesn't have same length")

    # Compare two sequence
    ref = ref.replace("E", "-")
    seq = seq.replace("E", "-")
    ops = []
    for i in range(len(ref)):
        if ref[i] == "-" and seq[i] == "-":
            continue
        if ref[i] == seq[i]:
            ops.append("M")
        elif ref[i] == "-":
            ops.append("I")
        elif seq[i] == "-":
            ops.append("D")
        else:
            ops.append("X")

    # aggregate op
    # ops: MMMIIDDD
    # -> cigar: M3I2D3
    op_type, count = "", 0
    cigar = []
    for op in ops:
        if op != op_type:
            if count:
                cigar.append((op_type, count))
            count = 0
        op_type = op
        count += 1
    if count:
        cigar.append((op_type, count))

    return cigar
