# pyHLAMSA

This is a small python utility that deal with MSA(multiple sequence alignment)

Implemented with some useful functions shown below.

Still in development.

## Features

### 1. Read sequences from database in one line

`HLAmsa` provided a simple way to read HLA data

It can automatically download and read the sequences

``` python
>>> from pyHLAMSA import HLAmsa

>>> hla = HLAmsa(["A", "B"], filetype="gen",
                 imgt_alignment_folder="alignment_v3470",
                 version="3470")

>>> print(hla.list_genes())
['A', 'B']

>>> print(hla["A"])
<A gen alleles=4101 block=5UTR(301) exon1(80) intron1(130) exon2(335) intron2(273) exon3(413) intron3(666) exon4(287) intron4(119) exon5(117) intron5(444) exon6(33) intron6(143) exon7(48) intron7(170) exon8(5) 3UTR(302)>
```

`KIRmsa` can also read sequences of KIR

``` python
>>> from pyHLAMSA import KIRmsa

# If don't specific the genes, it will read all genes.
>>> kir = KIRmsa(ipd_folder="KIR_v2110", version="2110")

>>> print(kir.list_genes())
['KIR2DL1', 'KIR2DL2', 'KIR2DL3', 'KIR2DL4', 'KIR2DL5', 'KIR2DP1', 'KIR2DS1', 'KIR2DS2', 'KIR2DS3', 'KIR2DS4', 'KIR2DS5', 'KIR3DL1', 'KIR3DL2', 'KIR3DL3', 'KIR3DP1', 'KIR3DS1']

>>> print(kir["KIR2DL1"])
<KIR2DL1 gen alleles=173 block=5UTR(268) exon1(34) intron1(964) exon2(36) intron2(728) exon3(282) intron3(1441) exon4(300) intron4(1534) exon5(294) intron5(3157) exon6(51) intron6(4270) exon7(102) intron7(462) exon8(53) intron8(98) exon9(177) 3UTR(510)>
```

### 2. Merge gene and nuc MSA is quiet simple

This main features give us a chance to use genomic MSA and nucleotide MSA at the same time.

The nucleotide MSA is a exon-only sequence, thus we fill the intron with `E` after merged.

``` python
>>> hla = HLAmsa(["A"], filetype=["gen", "nuc"],
                 imgt_alignment_folder="alignment_v3470")
>>> print(hla["A"])
<A gen alleles=7349 block=5UTR(301) exon1(80) intron1(130) exon2(351) intron2(273) exon3(436) intron3(666) exon4(361) intron4(119) exon5(117) intron5(444) exon6(33) intron6(143) exon7(48) intron7(170) exon8(5) 3UTR(302)>

# or manually
>>> a_gen = HLAmsa("A", filetype="gen", 
>>>                imgt_alignment_folder="alignment_v3470")["A"]

>>> print(a_gen)
<A gen alleles=4101 block=5UTR(301) exon1(80) intron1(130) exon2(335) intron2(273) exon3(413) intron3(666) exon4(287) intron4(119) exon5(117) intron5(444) exon6(33) intron6(143) exon7(48) intron7(170) exon8(5) 3UTR(302)>

>>> a_nuc = HLAmsa("A", filetype="nuc", 
>>>                imgt_alignment_folder="alignment_v3470")["A"]
>>> print(a_nuc)
<A nuc alleles=7353 block=exon1(80) exon2(351) exon3(436) exon4(361) exon5(117) exon6(33) exon7(48) exon8(5)>

>>> a_gen = a_gen.remove('A*03:437Q')
>>> print(a_gen.merge_exon(a_nuc))
<A gen alleles=7349 block=5UTR(301) exon1(80) intron1(130) exon2(351) intron2(273) exon3(436) intron3(666) exon4(361) intron4(119) exon5(117) intron5(444) exon6(33) intron6(143) exon7(48) intron7(170) exon8(5) 3UTR(302)>
```

### 3. Block-based selection: You can select ANY intron/exon.

``` python
# select exon2 and exon3
>>> exon23 = a_gen.select_exon([2,3])  # 1-base
>>> print(exon23)
<A nuc alleles=4100 block=exon2(335) exon3(413)>


# select exon2 + intron2 + exon3
>>> e2i2e3 = a_gen.select_block([3,4,5])  # 0-base
>>> print(e2i2e3)
<A  alleles=4100 block=exon2(335) intron2(273) exon3(413)>
```

### 4. Easy to compare alleles: select and print

``` python
# select first 10 alleles
>> exon23_10 = exon23.select_allele(exon23.get_sequence_names()[:10])
>>> print(exon23_10)
<A nuc alleles=10 block=exon2(335) exon3(413)>

# print it
>>> print(exon23_10.format_alignment())
                  812                                                      1136
                    |                                                         |
 A*01:01:01:01      ACCGAGCGAA CCTGGGGACC CTGCGCGGCT ACTACAACCA GAGCGAGGAC G| GTTCTCACA CC-ATCCAGA TAATGTATGG CTGCGACG-- ----------
 A*01:01:01:02N     ACCGAGCGAA CCTGGGGACC CTGCGCGGCT ACTACAACCA GAGCGAGGAC G| GTTCTCACA CC-ATCCAGA TAATGTATGG CTGCGACG-- ----------
 A*01:01:01:03      ACCGAGCGAA CCTGGGGACC CTGCGCGGCT ACTACAACCA GAGCGAGGAC G| GTTCTCACA CC-ATCCAGA TAATGTATGG CTGCGACG-- ----------
 A*01:01:01:04      ACCGAGCGAA CCTGGGGACC CTGCGCGGCT ACTACAACCA GAGCGAGGAC G| GTTCTCACA CC-ATCCAGA TAATGTATGG CTGCGACG-- ----------
 A*01:01:01:05      ACCGAGCGAA CCTGGGGACC CTGCGCGGCT ACTACAACCA GAGCGAGGAC G| GTTCTCACA CC-ATCCAGA TAATGTATGG CTGCGACG-- ----------
 A*01:01:01:06      ACCGAGCGAA CCTGGGGACC CTGCGCGGCT ACTACAACCA GAGCGAGGAC G| GTTCTCACA CC-ATCCAGA TAATGTATGG CTGCGACG-- ----------
 A*01:01:01:07      ACCGAGCGAA CCTGGGGACC CTGCGCGGCT ACTACAACCA GAGCGAGGAC G| GTTCTCACA CC-ATCCAGA TAATGTATGG CTGCGACG-- ----------
 A*01:01:01:08      ACCGAGCGAA CCTGGGGACC CTGCGCGGCT ACTACAACCA GAGCGAGGAC G| GTTCTCACA CC-ATCCAGA TAATGTATGG CTGCGACG-- ----------
 A*01:01:01:09      ACCGAGCGAA CCTGGGGACC CTGCGCGGCT ACTACAACCA GAGCGAGGAC G| GTTCTCACA CC-ATCCAGA TAATGTATGG CTGCGACG-- ----------
 A*01:01:01:10      ACCGAGCGAA CCTGGGGACC CTGCGCGGCT ACTACAACCA GAGCGAGGAC G| GTTCTCACA CC-ATCCAGA TAATGTATGG CTGCGACG-- ----------

# using regex to select
# "|" indicate the border of block, in this case, it's the border of exon2 and exon3
>>> exon23_1field = exon23.select_allele(r"A\*.*:01:01:01$")
>>> print(exon23_1field.format_alignment_diff())
                  812                                                      1136
                    |                                                         |
 A*01:01:01:01      ACCGAGCGAA CCTGGGGACC CTGCGCGGCT ACTACAACCA GAGCGAGGAC G| GTTCTCACA CC.ATCCAGA TAATGTATGG CTGCGACG.. ..........
 A*02:01:01:01      ------T-G- ---------- ---------- ---------- --------C- -| --------- --.G------ GG-------- --------.. ..........
 A*03:01:01:01      ------T-G- ---------- ---------- ---------- --------C- -| --------- --.------- ---------- --------.. ..........
 A*11:01:01:01      ------T-G- ---------- ---------- ---------- ---------- -| --------- --.------- ---------- --------.. ..........
 A*23:01:01:01      ------A--- ----C---T- GC--T-C--- ---------- --------C- -| --------- --.C------ -G----T--- --------.. ..........
 A*25:01:01:01      ------A--G ----C---T- GC--T-C--- ---------- ---------- -| --------- --.------- GG-------- --------.. ..........
 A*26:01:01:01      ---------- ---------- ---------- ---------- ---------- -| --------- --.------- GG-------- --------.. ..........
 A*29:01:01:01      ---------- ---------- ---------- ---------- --------C- -| --------- --.------- -G-------- ----C---.. ..........
 A*30:01:01:01      ------T-G- ---------- ---------- ---------- --------C- -| --------- --.------- ---------- --------.. ..........
 A*32:01:01:01      ------A--G ----C---T- GC--T-C--- ---------- --------C- -| --------- --.------- -G-------- --------.. ..........
 A*33:01:01:01      ------T-G- ---------- ---------- ---------- --------C- -| --------- --.------- -G-------- --------.. ..........
 A*34:01:01:01      ------T-G- ---------- ---------- ---------- ---------- -| --------- --.------- GG-------- --------.. ..........
 A*36:01:01:01      ---------- ---------- ---------- ---------- ---------- -| --------- --.------- ---------- --------.. ..........
 A*66:01:01:01      ------T-G- ---------- ---------- ---------- ---------- -| --------- --.------- GG-------- --------.. ..........
 A*68:01:01:01      ------T-G- ---------- ---------- ---------- --------C- -| --------- --.------- -G-------- --------.. ..........
 A*69:01:01:01      ------T-G- ---------- ---------- ---------- --------C- -| --------- --.G------ GG-------- --------.. ..........
 A*74:01:01:01      ------T-G- ---------- ---------- ---------- --------C- -| --------- --.------- -G-------- --------.. ..........
 A*80:01:01:01      ---------- ---------- ---------- ---------- ---------- -| --------- --.------- ---------- --------.. ..........

# show only variantiation
>>> print(exon23_1field.format_variantion_base())
Total variantion: 71
                         536 537  541          565 567 570          593          640          654 658          684          721      
                           |   |    |            |   |   |            |            |            |   |            |            |      
 A*01:01:01:01      | ATTTCT   TCAC ATCCG| CCGGC CG  CGG GGA.G | ATCGCCGTGG | .G.ACACG.C | CGTGCGGTTC GACA| AGA..A GATG| CGGGCG CCGT|
 A*02:01:01:01      | ------   ---- -----| ----- --  --- ---.- | -----A---- | .-.-----.- | ---------- ----| ---..G ----| ------ ----|
 A*03:01:01:01      | ------   ---- -----| ----- --  --- ---.- | ---------- | .-.-----.- | ---------- ----| ---..G ----| ------ ----|
 A*11:01:01:01      | ------   A--- C----| ----- --  --- ---.- | ---------- | .-.-----.- | ---------- ----| ---..G ----| ------ ----|
 A*23:01:01:01      | ------   C--- -----| ----- --  --- ---.- | ---------- | .-.-----.- | ---------- ----| ---..G ----| ------ ----|
 A*25:01:01:01      | ------   A--- C----| ----- --  --- ---.- | ---------- | .-.-----.- | ---------- ----| ---..G ----| ------ ----|
 A*26:01:01:01      | ------   A--- C----| ----- --  --- ---.- | ---------- | .-.-----.- | ---------- ----| ---..G ----| ------ ----|
 A*29:01:01:01      | -----A   C--- -----| ----- --  --- ---.- | ---------- | .-.-----.- | ---------T ----| ---..G ----| -----A ----|
 A*30:01:01:01      | ------   C--- -----| ----- A-  T-- A--.- | -----A---- | .-.-----.- | ---------- ----| ---..G ----| ------ ----|
 A*32:01:01:01      | ------   ---- -----| ----- --  --- ---.- | ---------- | .-.-----.- | ---------T ----| ---..G ----| ------ ----|
 A*33:01:01:01      | -----A   C--- -----| ----- --  --- ---.- | ---------- | .-.-----.- | ---------- ----| ---..G ----| ------ ----|
 A*34:01:01:01      | ------   A--- C----| ----- --  --- ---.- | ---------- | .-.-----.- | ---------- ----| ---..G ----| ------ ----|
 A*36:01:01:01      | ------   ---- -----| ----- --  --- ---.- | ---------- | .-.-----.- | ---------- ----| ---..- ----| ------ ----|
 A*66:01:01:01      | ------   A--- C----| ----- --  --- ---.- | ---------- | .-.-----.- | ---------- ----| ---..G ----| ------ ----|
 A*68:01:01:01      | ------   A--- C----| ----- --  --- ---.- | ---------- | .-.-----.- | ---------- ----| ---..G ----| ------ ----|
 A*69:01:01:01      | ------   A--- C----| ----- --  --- ---.- | ---------- | .-.-----.- | ---------- ----| ---..G ----| ------ ----|
 A*74:01:01:01      | ------   ---- -----| ----- --  --- ---.- | ---------- | .-.-----.- | ---------T ----| ---..G ----| ------ ----|
 A*80:01:01:01      | ------   ---- -----| ----- --  --- ---.- | -----A---- | .-.--T--.- | -----A---- ----| ---..G ----| ------ ----|
```

### 5. Some useful operation

``` python
# Calculate variant frequency of ( A T C G - ) per base
>>> print(exon23.calculate_frequency()[:10])
[[2, 3, 0, 4095, 0], [2, 1, 4097, 0, 0], [0, 4098, 2, 0, 0], [0, 3, 4095, 2, 0], [0, 911, 3188, 0, 1], [0, 1, 4097, 1, 1], [0, 0, 1, 0, 4099], [4097, 0, 0, 2, 1], [4, 3, 4090, 3, 0], [1, 4097, 1, 1, 0]]

# Get consensus(The largest one in frequency) among the msa
>>> exon23.get_consensus(include_gap=True)
'GCTCCC-ACTCCATGAGGTATTTCTTCACATCCGTGTCCCGGCCCGGCCGCGGGGA----GCCCCGCTTCATCGCCGTGGGC-----------------------TACGTGGACG-ACACG-CAGTTCGTGCGGTTCGACAGCGACGCCGCGAGCCAGAGGATGGAGCCG--------------------CGGGCGCCGTGGATA-GAGCAGGAGGGGCCGGAGTATTGGGACCAGGAGACACGGA-------------A-TGTGAAGGCCCACTCACAGACTGACCGAGTGGACCTGGGGACCCTGCGCGGCTACTACAACCAGAGCGAGGCCGGTTCTCACACC-ATCCAGATGATGTATGGCTGCGACG--------------TGGGG-TCGGACGGGCGCTTCCTCCGCGGGTACCAGCA---GGACGCCTACGACGGCAAGGATTAC---ATCGCCCTGAAC------------------------GAGGACCTGCGCTCTTGGACCGCGGCGGAC--------ATGGCGGCTCAGATCACCAAGCGC-AAGT----GGGAGG--CGGCCC-ATGT------------------------------------------GGCGG-AGCAGTTGAGAGCCTACCTGGAGGGCACG--------TGCGTG----GAGTGGCTCCG--CAGATA-CCTGGAGAACGGGAAGGAGACGCTGCAGC-----------------GCACGG'

# Add sequence into MSA
# include_gap=False in get_consensus will ignore the frequency of gap. i.e. choose one of ATCG
>>> a_merged = hla["A"]
>>> consensus_seq = a_merged.get_consensus(include_gap=False)
>>> a_merged.append("A*consensus", consensus_seq)
>>> a_merged.fill_imcomplete("A*consensus")

# Shrink: remove gap if all bases in the column are gap
>>> print(exon23_10.shrink().format_alignment_diff())
 gDNA               200
                    |
 A*01:01:01:01      AAGGCCCACT CACAGACTGA CCGAGCGAAC CTGGGGACCC TGCGCGGCTA CTACAACCAG AGCGAGGACG| GTTCTCACAC CATCCAGATA ATGTATGGCT
 A*01:01:01:02N     ---------- ---------- ---------- ---------- ---------- ---------- ----------| ---------- ---------- ----------
 A*01:01:01:03      ---------- ---------- ---------- ---------- ---------- ---------- ----------| ---------- ---------- ----------
 A*01:01:01:04      ---------- ---------- ---------- ---------- ---------- ---------- ----------| ---------- ---------- ----------
 A*01:01:01:05      ---------- ---------- ---------- ---------- ---------- ---------- ----------| ---------- ---------- ----------
 A*01:01:01:06      ---------- ---------- ---------- ---------- ---------- ---------- ----------| ---------- ---------- ----------
 A*01:01:01:07      ---------- ---------- ---------- ---------- ---------- ---------- ----------| ---------- ---------- ----------
 A*01:01:01:08      ---------- ---------- ---------- ---------- ---------- ---------- ----------| ---------- ---------- ----------
 A*01:01:01:09      ---------- ---------- ---------- ---------- ---------- ---------- ----------| ---------- ---------- ----------
 A*01:01:01:10      ---------- ---------- ---------- ---------- ---------- ---------- ----------| ---------- ---------- ----------

# select specific column 
>>> a_gen[12:100]
<A  alleles=4100 block=(88)>

>>> print(a_gen[[12,100]].format_alignment())
 gDNA               0
                    |
 A*01:01:01:01      GG
 A*01:01:01:02N     -G
 A*01:01:01:03      GG
 A*01:01:01:04      --
 A*01:01:01:05      --
 A*01:01:01:06      --
 A*01:01:01:07      --
 A*01:01:01:08      --
 A*01:01:01:09      --
 A*01:01:01:10      --
 A*01:01:01:11      GG
 A*01:01:01:12      GG

# concat
>>> print(a_gen.select_exon([2]) + a_gen.select_exon([3]))
<A nuc alleles=4100 block=exon2(335) exon3(413)>
```


### 6. Write and export the MSA

**Causion: If you run `merge_exon`, or the msa is from `filetype=['gen', 'nuc']`,
You should fill the `E` BEFORE save it.**

You can fill it by consensus_seq shown before.


* IMGT alignment format

``` python
print(a_gen.format_alignment_diff(), file=open("filename.txt", "w"))
```

* MultipleSeqAlignment

Transfer to [MultipleSeqAlignment](https://biopython.org/docs/1.75/api/Bio.Align.html#Bio.Align.MultipleSeqAlignment)

``` python
>>> print(a_gen.to_MultipleSeqAlignment())
Alignment with 4100 rows and 3866 columns
CAGGAGCAGAGGGGTCAGGGCGAAGTCCCAGGGCCCCAGGCGTG...AAA A*01:01:01:01
--------------------------------------------...--- A*01:01:01:02N
CAGGAGCAGAGGGGTCAGGGCGAAGTCCCAGGGCCCCAGGCGTG...AAA A*01:01:01:03
--------------------------------------------...--- A*01:01:01:04
--------------------------------------------...AAA A*01:01:01:05
--------------------------------------------...--- A*01:01:01:06
--------------------------------------------...AAA A*01:01:01:07
--------------------------------------------...--- A*01:01:01:08
--------------------------------------------...AAA A*01:01:01:09
--------------------------------------------...--- A*01:01:01:10
CAGGAGCAGAGGGGTCAGGGCGAAGTCCCAGGGCCCCAGGCGTG...--- A*01:01:01:11
CAGGAGCAGAGGGGTCAGGGCGAAGTCCCAGGGCCCCAGGCGTG...AAA A*01:01:01:12
--------------------------------------------...AAA A*01:01:01:13
--------------------------------------------...AAA A*01:01:01:14
--------------------------------------------...--- A*01:01:01:15
CAGGAGCAGAGGGGTCAGGGCGAAGTCCCAGGGCCCCAGGCGTG...--- A*01:01:01:16
CAGGAGCAGAGGGGTCAGGGCGAAGTCCCAGGGCCCCAGGCGTG...--- A*01:01:01:17
CAGGAGCAGAGGGGTCAGGGCGAAGTCCCAGGGCCCCAGGCGTG...AAA A*01:01:01:18
...
--------------------------------------------...--- A*80:07
```

* list of SeqRecord
``` python
# Save as MSA
SeqIO.write(a_gen.to_fasta(gap=True), "filename.msa.fa", "fasta")
# Save as no-gapped sequences
SeqIO.write(a_gen.to_fasta(gap=False), "filename.fa", "fasta")
```

* bam
``` python
a_gen.save_bam("filename.bam")
```

* gff
``` python
a_gen.save_gff("filename.gff")
```

After save the MSA as bam and gff, you can show the alignments on IGV
![msa_igv_example](https://raw.githubusercontent.com/linnil1/pyHLAMSA/main/HLA_msa.png)

* save/load
I use json and fasta to save our model
``` python
from pyHLAMSA import Genemsa
a_gen.save_msa("a_gen.fa", "a_gen.json")
a_gen = Genemsa.load_msa("a_gen.fa", "a_gen.json")
```


## TODO
* [x] Testing
  * [x] Main function
  * [x] exon-only sequence handling
  * [x] Reading from file
* [x] Some useful function: `copy`, `remove`, `get_sequence_names`, `__len__`, `size`, `split`, `concat`
* [x] Cannot handle pseudo exon
* [x] merge blocks and labels
* [x] Fix KIR when merge gen and nuc
* [x] Use index to trace orignal position
* [ ] Download latest version of IMGT or IPD
* [ ] CYP


## Requirement
* python3.8
* biopython
* pysam
* wget
* git


## Installation
```
git clone https://github.com/linnil1/pyHLAMSA
pip3 install -e pyHLAMSA
```


## Setup Document
```
pip3 install mkdocs mkdocs-material mkdocstrings
mkdocs serve
```


## Test
```
pip3 install -e .
pip3 install pytest pycodestyle
pytest
pycodestyle pyHLAMSA/*.py
```


## Some QA
> Why not inherit Bio.AlignIO.MultipleSeqAlignment?

The class does't support lot of functions than I expected.

And it's tidious to overwrite most of the functions to maintain our blocks information

> Why not use numpy to save the sequence?

Performance issue is not my bottle-neck yet.


## Citation
* IGV
    James T. Robinson, Helga Thorvaldsdóttir, Wendy Winckler, Mitchell Guttman, Eric S. Lander, Gad Getz, Jill P. Mesirov. Integrative Genomics Viewer. Nature Biotechnology 29, 24–26 (2011)
* IMGT
    Robinson J, Barker DJ, Georgiou X, Cooper MA, Flicek P, Marsh SGE
    IPD-IMGT/HLA Database
    Nucleic Acids Research (2020) 48:D948-55
* This github

## Document
See [https://linnil1.github.io/pyHLAMSA](See https://linnil1.github.io/pyHLAMSA)

::: pyHLAMSA
