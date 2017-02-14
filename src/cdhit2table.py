## !/usr/bin/env python
"""
Consider that last ID is separated as samples_id
"""
from Bio import SeqIO
import pandas as pd
import numpy as np
import click
from os import path


@click.command()
@click.option("-cdhit_out", help="CD-HIT ouput as fasta file", type=str,
              default="./test.fa", show_default=True)
@click.option("-out", help="Output file name", type=str, default="../test.tab",
              show_default=True)
@click.option("-sorttab",
              help="Short table based on sample and sequence count",
              type=bool, default=True, show_default=True)
def run(cdhit_out, out, sorttab):
    """Convert cdhit output to table."""
    if not cdhit_out:
        click.echo("Imputfile is not given. Exiting ...")
        exit(1)
    print(cdhit_out)
    if not path.isfile(cdhit_out):
        click.echo("Given input file doesn't exist. Exiting ....")
        exit(1)
    try:
        SeqIO.parse(cdhit_out, "fasta")
        table = {'cluster':[]}
        seq_sizes = {}
        current_cluster = 0
        seq_count = {}
        samp_count = {}
        with open(cdhit_out) as fin:
            for line in fin:
                if line[0] == '>':
                    cluster = line[1:-1].replace(" ", "_")
                    current_cluster += 1
                    seq_sizes[cluster] = []
                    table["cluster"].append(cluster)
                    continue
                data = line.split()
                sequenceid = data[2][1:-3]
                seqsize = int(data[1][:-3])
                seq_sizes[cluster].append(seq_sizes)
                samp = '_'.join(sequenceid.split("_")[:-1])
                if samp not in table:
                    # TODO : A lot 
                    table[samp] = ['*'] * (current_cluster - 1)
                    table[samp].append(sequenceid)
                    # include the possibilities of multiple sequences
                print(data)
    except IOError:
        click.echo("Given file is not in Fasta format. Exiting ...")
        exit(1)

    pass


if __name__ == '__main__':
    run()
