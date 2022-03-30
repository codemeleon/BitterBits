#!/usr/bin/env python

import click
import pandas as pd
from Bio import SeqIO


@click.command()
@click.option(
    "--aln",
    help="Fasta alignment file",
    type=click.File(),
    default=None,
    required=True,
    show_default=True,
)
@click.option(
    "--tab",
    help="Table with ids and year columns",
    type=click.File(),
    default=None,
    required=True,
    show_default=True,
)
@click.option(
    "--out_file",
    help="Alignment file containing year informaion in the sequence ids",
    type=click.File(mode="w"),
    default=None,
    required=True,
    show_default=True,
)
def run(aln, tab, out_file):
    """Add year information to a fasta alignment file."""
    year = pd.read_csv(tab)
    for rec in SeqIO.parse(aln, "fasta"):
        year_id = year.loc[year["id"] == rec.id, "year"].values[0]
        rec.id = rec.id + "_" + str(year_id)
        out_file.write(">" + rec.id + "\n" + str(rec.seq) + "\n")


if __name__ == "__main__":
    run()
