#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File              : gene_aln_blat.py
# Date              : 13.08.2021
# Last Modified Date: 13.08.2021
# -*- coding: utf-8 -*-
# File              : gene_aln_blat.py
# Date              : 13.08.2021
# Last Modified Date: 13.08.2021

from Bio import SeqIO
import click
import pandas as pd
from os import system


def run(ref_seq, assemblies):
    """Generate alignment of given gene along with reference and sequences from
    different samples

    :ref_seq: file containing ref gene
    :assemblies: denovo genome assembly by mapping to reference genome.
    :returns: None

    """

    if not ref_seq:
        exit("ref_seq not given")
    if not assemblies:
        exit("assemblies file not given")
    if path.exists(ref_seq):
        exit()

    refs = {}
    for rec in SeqIO.parse(ref_seq, "fasta"):
        refs[rec.id] = rec.seq
    assembly = {}
    for rec in SeqIO.parse(assemblies, "fasta"):
        assemblies[rec.id] = rec.seq

        pass


if __name__ == "__main__":
    run()
