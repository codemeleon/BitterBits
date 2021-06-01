#!/usr/bin/env python

# Code to count bases
import base64

from sys import argv

from Bio import SeqIO

bases = []

for rec in SeqIO.parse(argv[1], "fasta"):
    bases.append(str(rec.seq[int(argv[2])]).upper())


bases = {x:bases.count(x) for x in bases}
print(bases)
