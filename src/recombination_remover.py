#!/usr/bin/env python

import click
import pandas as pd
import numpy as np
from Bio import SeqIO


def merged_range(intervals):
    sorted_by_lower_bound = sorted(intervals, key=lambda tup: tup[0])
    merged = []

    for higher in sorted_by_lower_bound:
        if not merged:
            merged.append(higher)
        else:
            lower = merged[-1]
            # test for intersection between lower and higher:
            # we know via sorting that lower[0] <= higher[0]
            if higher[0] <= lower[1]:
                upper_bound = max(lower[1], higher[1])
                # replace by merged interval
                merged[-1] = (lower[0], upper_bound)
            else:
                merged.append(higher)
    return merged


@click.command()
@click.option(
    "--gff",
    help="Gubbin reconbination prediction gff file",
    type=click.File(),
    default=None,
    required=True,
    show_default=True,
)
@click.option(
    "--aln",
    help="Alignement file",
    type=click.File(),
    default=None,
    required=True,
    show_default=True,
)
@click.option(
    "--out_file",
    help="Output alignment file after removal of recombinations spot",
    type=click.File(mode="w"),
    default=None,
    required=True,
    show_default=True,
)
def run(aln, gff, out_file):
    """"""
    df = pd.read_table(gff, comment="#", header=None, usecols=[3, 4])
    ranges = merged_range(df.drop_duplicates().values.tolist())
    ranges.reverse()
    # with open(out_file, "w") as fout:
    for rec in SeqIO.parse(aln, "fasta"):
        tseq = str(rec.seq)
        for start, end in ranges:
            tseq = tseq[:start] + tseq[end:]
        out_file.write(">{}\n{}\n".format(rec.id, tseq))


if __name__ == "__main__":
    run()
