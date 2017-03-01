##!/usr/bin/env python

"""Generate fasta files for clusters."""

from Bio import SeqIO
import pandas as pd
from glob import glob
import click
from os import path


@click.command()
@click.option("-faaf", help="Fast sequence folder",
              type=str, default=None, show_default=True)
@click.option("-tabfile", help="Cluster Table file",
              type=str, default=None, show_default=True)
@click.option("-outf", help="Ouput folder for cluster fastas",
              type=str, default=None, show_default=True)
@click.option("-minseq",
              help="Minimumn number of sequences present in cluster",
              type=int, default=2, show_default=True)
def run(faaf, tabfile, outf, minseq):
    """"Generate Cluster fasta files base on given table."""

    assert path.isdir(faaf), "Input sequence folder not give. Exiting ..."
    assert path.isfile(tabfile), "Input table file not give. Exiting ..."
    assert path.isfile(tabfile), "Output cluster sequencefolder not give.\
        Exiting ..."
    table = pd.read_table(tabfile, index_col="cluster")
    table = table[table["seq_count"] >= minseq]
    for col in ["samp_count", "seq_count", "min", "median",
                "mean", "std", "max"]:
        del table[col]
    sequences = {}
    for fl in glob("%s/*" % faaf):
        sequences.update(SeqIO.to_dict(SeqIO.parse(fl, "fasta")))
    for cluster, row in table.iterrows():
        outfile = open("%s/%s.fasta" % (outf, cluster), "w")
        for seq_ids_len in row:
            if seq_ids_len == "*":
                continue
            seq_ids_len = seq_ids_len.split(";")
            for seq_id_len in seq_ids_len:
                seq_id = seq_id_len.split(":")[0]
                outfile.write(">%s\n%s\n" % (seq_id, sequences[seq_id].seq))
        outfile.close()


if __name__ == '__main__':
    run()
