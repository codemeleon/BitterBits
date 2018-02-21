from os import path

import click
import pandas as pd
import regex as re
from Bio import SeqIO


def longest_seq(seq):
    startP = re.compile('ATG')
    longest = (0, )
    for m in startP.finditer(str(seq), overlapped=True):
        if len(seq[m.start():].translate(to_stop=True)) > longest[0]:
            pro = seq[m.start():].translate(to_stop=True)
            longest = (len(pro), m.start(), str(pro),
                       seq[m.start():m.start() + len(pro) * 3 + 3])
    return longest[-1]
    # Add other details


@click.command()
@click.option(
    "-fi",
    help="Fasta file in",
    type=str,
    default="/home/devil/Documents/Tools/BitterBits/data/alntrimmer/FluB_iva_HA.fa",
    show_default=True)
@click.option(
    "-fo", help="Fasta file out", type=str, default=None, show_default=True)
def run(fi, fo):
    """This script is to extract the sequences for the logest frame on
    the forward strand of the sequences"""
    if not fi:
        exit("Input file not given. Exiting . . . .")
    if not path.isfile(fi):
        exit("Give path either doesn't exist or it is not a file")
    try:
        fasta = SeqIO.parse(fi, "fasta")
    except ValueError:
        exit("Given file is not in fasta format")

    sequences = {}
    for rec in fasta:
        sequences[rec.id] = longest_seq(rec.seq)
    if fo:
        with open(fo, "w") as fout:
            for k in sequences:
                fout.write(">%s\n%s\n" % (k, sequences[k]))
    else:
        for k in sequences:
            print(">%s\n%s\n" % (k, sequences[k]))


if __name__ == '__main__':
    run()
