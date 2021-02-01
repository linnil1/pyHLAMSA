# pyHLAMSA

This is a small python utility to parse alignments in IMGT.

It implemented with some useful functions shown below.

Still in development.

## Feature

### Download IMGT and read
It can be automatically download and parse IMGAHLA gen and nuc.

``` python
# A simple interface read A B DPA1 allele
# If the txt file not exist in imgt_folder, it will download itself
hla = HLAmsa(["A", "B", "DPA1"], filetype=["gen", "nuc"], imgt_folder="alignments", version="3430")
a = hla.genes["A"]
print(a)


# Genemsa is a main class to store alignments and implemented some important function
a = Genemsa("A", "gen")
a.read_file("alignments/A_gen.txt")
print(a)
```

### Merging
Merge gen and nuc alignments is easy

``` python
# If gen and nuc are both in filetype, it will merge
hla = HLAmsa(["A"], filetype=["gen", "nuc"])
# It's ok to read one of them without mergin
hla1 = HLAmsa(["A"], filetype=["gen"])
hla2 = HLAmsa(["A"], filetype=["nuc"])


# Or merge by yourself
a_gen = Genemsa("A", "gen")
a_gen.read_file("alignments/A_gen.txt")
a_nuc = Genemsa("A", "nuc")
a_nuc.read_file("alignments/A_nuc.txt")
a_gen = a_gen.merge_exon(a_nuc)
```

### Consensus
Calculate variant frequency and consensus

``` python
print(a.calculate_frequency())
consensus_seq = a.get_consensus(include_gap=False)

# fill the exon-only allele by consensus
a.add("A*consensus", consensus_seq)
a.fill_imcomplete("A*consensus")
```

### Selection
Select and insert sequences

``` python
# use regex to choose what allele you want
a_sub = a.select_allele(r"DPA1\*.*:01:01:01$")

# insert one sequence
a_sub.add("A*consensus", a.get_consensus(include_gap=False))

# extract the alignments from 50 to 199 bp
a_sub = a_sub[50:200]

# extract exon2 and exon3
a_sub.select_exon([2,3])

# extract exon2, intron2, exon3
a_sub.select_chunk([3,4,5])
```

### Output
Write to msf or alignment or other format

``` python
# print object
print(a_sub)

# print one sequence
a_sub.alleles["query"]

# print raw alignments
print(a_sub.select_exon().format_alignment())

# print diff alignment(Look like xx_gen.txt)
print(a_sub.select_exon([]).format_alignment_diff("query"))

# Convert object to MultipleSeqAlignment(Bio.Align)
# Thus, you can save in any format
print(a_sub.to_MultipleSeqAlignment())

# save to fasta(no gap)
SeqIO.write(a_sub.to_fasta(gap=False), "tmp.fa", "fasta")

# save to bam file
a_sub.save_bam("tmp.bam", ref_allele="A*consensus")

# save to gff3 (This file can show where exons are in IGV)
a_sub.save_gff("tmp.gff")
```

### Example

see main in `HLAmsa.py`

```
 A*consensus        ATGGCCGTCA TGGCGCCCCG AACCCTCGTC CTGCTACTCT CGGGGGCCCT GGCCCTGGCC CTGACCCAGA CCTGGGCGGG| GCTCCCCACT CCATGAGGTA
 A*01:01:01:01      ---------- ---------- -------C-- ---------- ---------- ------**** **-------- ---------*| ------*--- ----------
 A*02:01:01:01      ---------- ---------- ---------- ---------- -------T-- ------**** **-------- ---------*| ----T-*--- ----------
 A*03:01:01:01      ---------- ---------- -------C-- ---------- ---------- ------**** **-------- ---------*| ------*--- ----------
 A*11:01:01:01      ---------- ---------- -------C-- ---------- ---------- ------**** **-------- ---------*| ------*--- ----------
 A*23:01:01:01      ---------- ---------- ---------- ---------- ---------- ------**** **-------- -------A-*| ------*--- ----------
 A*25:01:01:01      ---------- ---------- ---------- ---------- ---------- ------**** **-------- ---------*| ------*--- ----------
 A*26:01:01:01      ---------- ---------- ---------- ---------- ---------- ------**** **-------- ---------*| ------*--- ----------
 A*29:01:01:01      ---------- ---------- -------C-- ---------- T--------- ------**** **-------- ---------*| ------*--- ----------
 A*30:01:01:01      ---------- ---------- -------C-- ---------- ---------- ------**** **-------- ---------*| ------*--- ----------
 A*32:01:01:01      ---------- ---------- -------C-- ---------- T--------- ------**** **-------- ---------*| ------*--- ----------
 A*33:01:01:01      ---------- ---------- -------C-- ---------- T--------- ------**** **-------- ---------*| ------*--- ----------
 A*34:01:01:01      ------A--- ---------- ---------- ---------- ---------- ------**** **-------- ---------*| ------*--- ----------
 A*36:01:01:01      ---------- ---------- -------C-- ---------- ---------- ------**** **-------- ---------*| ------*--- ----------
 A*66:01:01:01      ---------- ---------- ---------- ---------- ---------- ------**** **-------- ---------*| ------*--- ----------
 A*68:01:01:01      ---------- ---------- ---------- ---------- ---------- ------**** **-------- ---------*| ------*--- ----------
 A*69:01:01:01      ---------- ---------- ---------- ---------- ---------- ------**** **-------- ---------*| ------*--- ----------
 A*74:01:01:01      ---------- ---------- -------C-- ---------- T--------- ------**** **-------- --A------*| ------*--- ----------
 A*80:01:01:01      ---------- --C------- -------C-- ---------- ---------- ------**** **-------- -------A-*| ------*--- ----------
```

You can show the alignments on IGV
![](https://raw.githubusercontent.com/linnil1/pyHLAMSA/main/HLA_msa.png)


## TODO
* [ ] Sanity check
* [ ] Cannot handle splice variant


## Requirement
* python3.8
* biopython
* pysam
* wget

## Citation
* IGV
    James T. Robinson, Helga Thorvaldsdóttir, Wendy Winckler, Mitchell Guttman, Eric S. Lander, Gad Getz, Jill P. Mesirov. Integrative Genomics Viewer. Nature Biotechnology 29, 24–26 (2011)
* IMGT
    Robinson J, Barker DJ, Georgiou X, Cooper MA, Flicek P, Marsh SGE
    IPD-IMGT/HLA Database
    Nucleic Acids Research (2020) 48:D948-55
* This github
