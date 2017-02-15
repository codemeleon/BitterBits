#!/usr/bin/env python

"""Consider that ID is separated as samples_id."""
from Bio import SeqIO
import pandas as pd
import numpy as np
import click


@click.command()
@click.option("-cdhit_out", help="CD-HIT ouput as fasta file", type=str,
              default=None, show_default=True)
@click.option("-out", help="Output file name", type=str, default=None,
              show_default=True)
@click.option("-sorttab",
              help="Short table based on sample and sequence count in\
              decreasing order", type=bool, default=False, show_default=True)
def run(cdhit_out, out, sorttab):
    """Convert cdhit output to table."""
    # Expecting no empty line.
    assert cdhit_out, "Imputfile is not given. Exiting ..."
    assert cdhit_out, "Given input file doesn't exist. Exiting ...."
    assert SeqIO.parse(cdhit_out, "fasta"), "CDHIT output is not in Fasta.\
    Exiting ..."

    table = {}
    seq_sizes = []
    samp_count = []
    seq_min, seq_max, seq_mean, seq_median, seq_sd = [], [], [], [], []
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
    table["samp_count"] = samp_count
    table["seq_count"] = seq_count
    table["min"] = seq_min
    table["max"] = seq_max
    table["sd"] = seq_sd
    table["median"] = seq_median
    table["mean"] = seq_mean
    table = pd.DataFrame.from_dict(table)
    table["cluster"] = ["Cluster_%d" % i for i in range(len(table))]
    additiona_columns = ["cluster", "samp_count", "seq_count", "min",
                         "mean", "sd", "median", "max"]
    rearranged_columns = additiona_columns + list(set(table.columns) -
                                                  set(additiona_columns))
    table = table[rearranged_columns]
    if sorttab:
        table = table.sort_values(["samp_count", "seq_count"],
                                  ascending=[False, False])
    assert not table.to_csv(out, index=False, sep="\t"), "Location is write\
    protected. Exiting ..."


if __name__ == '__main__':
    run()
