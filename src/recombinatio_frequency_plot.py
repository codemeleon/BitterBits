#!/usr/bin/env python

import click
from os import path
import pandas as pd
from pylab import *


def ranges(p):
    """List to range."""
    # Copied from strackoverflow
    q = sorted(p)
    i = 0
    for j in range(1,len(q)):
        if q[j] > 1+q[j-1]:
            yield (q[i],q[j-1])
            i = j
    yield (q[i], q[-1])


@click.command()
@click.option("-gff", help="Input gff, generate by gubbins",
              type=str, default=None, show_default=True)
@click.option("-ofig", help="Output figure file, file format should be "
              "added in file extension", type=str, default=None,
              show_default=True)
@click.option("-fh", help="Figure height in cm", type=int,
              default=12, show_default=True)
@click.option("-fw", help="Figure width in cm", type=int,
              default=16, show_default=True)
def run(gff, ofig, fh, fw):
    """Generate recombination barplot from given gff file"""
    if not gff:
        click.echo("GFF file not given. Exiting ....")
        exit(1)
    if not ofig:
        click.echo("Output figure file not given. Exiting ....")
        exit(1)
    if not path.isfile(gff):
        click.echo("Given input file doesn't exist or it is not a file. "
                   "Exiting .....")
        exit(1)
    gff_table = pd.read_table(gff, header=None, comment="#", usecols=[3, 4, 8])
    gff_table['int_count'] = gff_table[8].map(lambda x: x.count('internal_'))
    del gff_table[8]
    gff_table = gff_table.loc[gff_table['int_count']==1]
    gff_table = gff_table.sort_values([3,4])
    gff_table = gff_table.drop_duplicates()

    pos = []
    for _, row in gff_table.iterrows():
        pos += list(range(row[3], row[4]+1))
    pos_count = {}
    pos.sort()
    for k in pos:
        if k in pos_count:
            pos_count[k] += 1
        else:
            pos_count[k] = 1
    del pos
    pos_count = pd.DataFrame.from_dict({'pos':list(pos_count.keys()),
                                        'count':list(pos_count.values())})
    pos_count_sorted = pos_count.sort_values("pos")
    

    x_range = []
    y = []
    for j in set(pos_count_sorted["count"]):
        rngs = list(ranges(pos_count_sorted.loc[
            pos_count_sorted["count"]==j, "pos"]))
        x_range += rngs
        y += [j] * len(rngs)

    for i, k in enumerate(x_range):
        fill_between(k, 0, [y[i]]*2, color='k')
    savefig(ofig)
    clf()

if __name__ == '__main__':
    run()
