#!/usr/bin/env python
from os import path

import click
from Bio import SeqIO


@click.command()
@click.option(
    '-fq', help="Input fastq", type=str, default=None, show_default=True)
@click.option(
    '-fa', help="Output fastq", type=str, default=None, show_default=True)
def run(fq, fa):
    """Converts fastq to fasta."""
    if not path.exists(fq) or not path.isfile(fq):
        exit(
            "Given fastq file doesn't exist or it is not a file. Exiting . . ."
        )
    if fa:
        outputfile = open(fa, "w")
    for rec in SeqIO.parse(fq, "fastq"):
        if fa:
            outputfile.write(">%s\n%s\n" % (rec.id, rec.seq))
        else:
            print(">%s\n%s\n" % (rec.id, rec.seq))
    if fa:
        outputfile.close()


if __name__ == '__main__':
    run()
