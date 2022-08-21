#!/usr/bin/env python

from Bio import SeqIO
import click


@click.command()
@click.option(
    "-i", "infile",
    help="Input fasta file.",
    type=str,
    default=None,
    show_default=True)
@click.option(
    "-o", "outfile",
    help="Output fasta file",
    type=str,
    default=None,
    show_default=True)
def run(infile, outfile):
    """
    Generates reverse complement of sequences in the given file.
    """
    with open(outfile) as fout:
        for rec in SeqIO.parse(infile, "fasta"):
            fout.write(f">{rec.description}\n{rec.seq.reverse_complement()}\n")


if __name__ == "__main__":
    run()
