#!/usr/bin/env python

"""Splits fasta sequencs from a file."""

from Bio import SeqIO
from os import path, makedirs
import click


@click.command()
@click.option("-infasta", help="Input fasta file", type=str,
              default=None, show_default=True)
@click.option("-outfold", help="Output folder", type=str,
              default=None, show_default=True)
@click.option("-id_include", help="Including id in file", type=bool,
              default=True, show_default=True)
def run(infasta, outfold, id_include):
    """Separate fasta sequence in difererent files."""
    assert infasta, "No input file given. Exiting ...."
    assert outfold, "No output folder given. Exiting ...."
    assert path.isfile(infasta), "Given input file doesn't exist. Exiting ..."
    assert SeqIO.parse(infasta, "fasta"), "Input is not fasta. Exiting ...."
    assert makedirs(outfold, exist_ok=True), "Unable to create given folder.\
        Exiting ...."
    if id_include:
        for rec in SeqIO.parse(infasta, "fasta"):
            with open("%s/%s.fasta" % (outfold, rec.id), "w") as fout:
                fout.write(">%s\n%s\n" % (rec.id, rec.seq))
    else:
        for rec in SeqIO.parse(infasta, "fasta"):
            with open("%s/%s.fasta" % (outfold, rec.id), "w") as fout:
                fout.write("%s\n" % rec.seq)


if __name__ == '__main__':
    run()
