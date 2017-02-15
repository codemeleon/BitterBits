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
    # Expecting no empty line.
    if not cdhit_out:
        click.echo("Imputfile is not given. Exiting ...")
        exit(1)
    print(cdhit_out)
    if not path.isfile(cdhit_out):
        click.echo("Given input file doesn't exist. Exiting ....")
        exit(1)
    try:
        SeqIO.parse(cdhit_out, "fasta")
        table = {}
        seq_sizes = []
        samp_count = []
        seq_min, seq_max, seq_mean, seq_median, seq_sd = [[]] * 5
        seq_count = []
        local_sets = {}
        current_cluster = 0
        samp_list = []
        with open(cdhit_out) as fin:
            for line in fin:
                if line[0] == '>':
                    for k in local_sets:
                        if k not in table:
                            table[k] = ['*'] * (current_cluster - 1)

                        table[k].append(local_sets[k])
                        # print(table)
                    for k in set(table.keys()) - set(local_sets.keys()):
                        table[k].append('*')
                    if len(seq_sizes):
                        samp_count.append(len(set(samp_list)))
                        samp_list = []
                        seq_count.append(len(seq_sizes))
                        seq_min.append(min(seq_sizes))
                        seq_max.append(max(seq_sizes))
                        seq_mean.append(np.mean(seq_sizes))
                        seq_median.append(np.median(seq_sizes))
                        seq_sd.append(np.std(seq_sizes))

                    seq_sizes = []
                    local_sets = {}
                    current_cluster += 1
                    continue
                data = line.split()
                sequenceid = data[2][1:-3]
                seqsize = int(data[1][:-3])
                seq_sizes.append(seqsize)
                samp = '_'.join(sequenceid.split("_")[:-1])
                samp_list.append(samp)
                if samp in local_sets:
                    local_sets[samp] += ",%s:%d" % (sequenceid, seqsize)
                else:
                    local_sets[samp] = "%s:%d" % (sequenceid, seqsize)
                    # include the possibilities of multiple sequences
            for k in local_sets:
                if k not in table:
                    table[k] = ['*'] * (current_cluster - 1)
                table[k].append(local_sets[k])

            for k in set(table.keys()) - set(local_sets.keys()):
                table[k].append('*')
        if seq_sizes:
            samp_count.append(len(set(samp_list)))
            seq_count.append(len(seq_sizes))
            seq_min.append(min(seq_sizes))
            seq_max.append(max(seq_sizes))
            seq_mean.append(np.mean(seq_sizes))
            seq_median.append(np.median(seq_sizes))
            seq_sd.append(np.std(seq_sizes))
            # local_sets = {}
        # print(len(samp_count), len(seq_count))
        table["samp_count"] = samp_count
        table["seq_count"] = seq_count
        table["min"] = seq_min
        table["max"] = seq_max
        table["sd"] = seq_sd
        table["median"] = seq_median
        table["mean"] = seq_mean
        # for k in table:
        print(table['mean'])
        # pd.DataFrame.from_dict(table)
    except IOError:
        click.echo("Given file is not in Fasta format. Exiting ...")
        exit(1)


if __name__ == '__main__':
    run()
