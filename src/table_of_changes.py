#!/usr/bin/env python

"""Generates table of changes based on aligned sequences."""

import click
from Bio import SeqIO
from glob import glob
import pandas as pd
from collections import Counter
import operator
from os import path


@click.command()
@click.option("-alnfold", help="Alignemnt folder path", type=str,
              default=None, show_default=True)
@click.option("-outfile", help="Output table file", type=str,
              default=None, show_default=True)
@click.option("-maxmaj", help="Max fraction of major varient", type=float,
              default=0.9, show_default=True)
def run(alnfold, outfile, maxmaj):
    """Generate table."""
    changes_df = pd.DataFrame()
    # cluster_seq_df = pd.DataFrame()
    for fl in glob("%s/*" % alnfold):
        t_df = {}
        # t_seq_df = {}
        for rec in SeqIO.parse(fl, "fasta"):
            t_df["_".join(rec.id.split("_")[:-1])] = list(str(rec.seq))
            # TODO: Modify split part based on your need in future
            # Currently based on  sanger's structure
        t_df = pd.DataFrame.from_dict(t_df)
        # print(t_df)
        # exit(1)
        selected_pos = []
        for i, row in t_df.iterrows():
            counts = dict(Counter(row))
            max_val_key = max(counts.items(), key=operator.itemgetter(1))[0]
            if counts[max_val_key] <= maxmaj * len(counts):
                selected_pos.append(i)
        t_df = t_df.ix[selected_pos]
        cluster = path.split(fl)[1].split(".")[0]
        t_df["cluster"] = [cluster] * len(t_df)
        changes_df = changes_df.append(t_df)
    changes_df.to_csv(outfile, index=False, sep="\t")


if __name__ == '__main__':
    run()
