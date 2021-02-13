from pyHLAMSA import HLAmsa, HLAmsaEX, Genemsa
from Bio import SeqIO
import logging

# setup logger to stdout
logger = logging.getLogger("pyHLAMSA")
logger.setLevel(logging.DEBUG)
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)
logger.addHandler(ch)

# Basic operation: read, add
# hla = HLAmsa(["A"], filetype=["gen", "nuc"],
#              imgt_folder="IMGT/alignments", version="3430")
hla = HLAmsaEX(["A"], filetype=["gen", "nuc"],
               imgt_folder="IMGT", version="3430")
a = hla.genes["A"]
a.add("A*consensus", a.get_consensus(include_gap=False))
a.fill_imcomplete("A*consensus")

# select and get consensus
a_sub = a.select_allele(r"A\*.*:01:01:01$")
a_sub.extend(a.select_allele(r"A\*consensus$"))
a_sub.sort()
print(a_sub)

# output
print(a_sub.select_exon().format_alignment_diff("A*consensus"))
print(a_sub.to_MultipleSeqAlignment())
a_sub.save_bam("tmp.bam", "A*consensus")
SeqIO.write(a_sub.to_fasta(gap=False), "tmp.fa", "fasta")
a_sub.save_gff("tmp.gff", strand="+")

# align seq on the consensus and print
seq = list(a.get("A*01:01:87"))
seq[55] = seq[65] = seq[78] = seq[79] = seq[80] = "T"
seq = a_sub.align("".join(seq), target_allele="A*consensus")
a_sub.add("A*query", seq)
print(a_sub[30:180].shrink().format_alignment_diff("A*consensus"))

# reverse the sequences if you needed
a_rv =  a_sub.reverse_complement()
a_rv.select_exon().save_gff("tmp_rv.gff", strand="-")
