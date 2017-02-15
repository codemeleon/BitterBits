##!/usr/bin/env python

"""Generate fasta files for clusters."""

from Bio import SeqIO
import pandas as pd
frm glob import glob
import click

@click.command()
@click.option("-faaf", help="Fast sequence folder",
              type=str, default=None, show_default=True)
@click.option("-tabfile", help="Cluster Table file",
              type=str, deafult=None, show_default=True)
@click.option("-outf", help="Ouput folder for cluster fastas",
              type=str, default=None, show_deafult=True)
@click.option("-minseq", help="Minimumn number of sequences present in cluster",
              type=int, default=2, show_default=True)
def run(faaf, tabfile, outf, minseq):
    pass




if  __name__ == '__main__':
    run()
