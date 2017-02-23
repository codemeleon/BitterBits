#!/usr/bin/env python
"""Licence GPLv3."""
import pandas as pd
from Bio import SeqIO
from glob import glob
from os import path, makedirs
import numpy as np
import click


def paralogs_separator(seq_ids, sequences_df):
    """Split paralogs."""
    max_seq_count = 0
    max_sq_list = []
    for k in seq_ids:
        if max_seq_count < len(seq_ids[k]):
            max_sq_list = seq_ids[k]
            max_seq_count = len(seq_ids[k])

    new_groups = {}
    for k in max_sq_list:
        new_groups[k] = [k]
    columns_to_compare = set(sequences_df.columns) - set(max_sq_list)

    for col in columns_to_compare:
        group = ''
        tdiff = 10000000
        for k in max_sq_list:
            diff = np.sum(sequences_df[k] != sequences_df[col])
            if diff < tdiff:
                tdiff = diff
                group = k
        new_groups[group].append(col)
    to_return = []
    for k in new_groups:
        qlist = []
        qdict = {}
        for l in new_groups[k]:
            tid = l.split('___')[0]
            qlist.append(tid)
            if tid in qdict:
                qdict[tid].append(l)
            else:
                qdict[tid] = [l]
        if len(set(qlist)) == len(qlist):
            to_return.append(new_groups[k])
        else:
            to_return += paralogs_separator(qdict,
                                            sequences_df[new_groups[k]])

    return to_return


# def testrecurse(z, target):
#     if z >= target:
#         return []
#     return [[z]] + testrecurse(2 * z, target)
# testrecurse(1, 1000)
# get the list of most common sample and sequences
#


@click.command()
@click.option("--ofile_folder", help="Orginal Sequence file folder",
              type=str, default=None, show_default=True)
@click.option("--cluster_aln_folder", help="Cluster Alignment folder",
              type=str, default=None, show_default=True)
@click.option("--outfolder", help="Separated paralog folder",
              type=str, default=None, show_default=True)
def run(ofile_folder, cluster_aln_folder, outfolder):
    """Split paralogs in different groups."""
    assert ofile_folder, "Original File folder not given. Exiting ..."
    assert cluster_aln_folder, "Cluster Alignment File folder not given.\
        Exiting ..."
    assert outfolder, "Output File folder not given. Exiting ..."
    assert path.isdir(ofile_folder), "Original file folder doesn't exit.\
        Exiting ..."
    assert path.isdir(cluster_aln_folder), "Cluster Alignment folder doesn't\
        exit. Exiting ..."
    assert not makedirs(outfolder, exist_ok=True), "Unable to create output folder.\
        Exiting ....."

    seq_samp = {}
    for fl in glob("%s/*" % ofile_folder):
        samp = path.split(fl)[1].split('.')[0]
        sequences = SeqIO.to_dict(SeqIO.parse(fl, "fasta"))
        for k in sequences.keys():
            seq_samp[k] = samp

    for fl in glob("%s/*" % cluster_aln_folder):
        file_id = path.split(fl)[1].split('.')[0]
        sequences_df = {}
        sequences = {}
        seq_ids = {}
        samples = []
        for rec in SeqIO.parse(fl, "fasta"):
            sequences_df["%s___%s" %
                         (seq_samp[rec.id], rec.id)
                         ] = list(str(rec.seq))

            sequences["%s___%s" %
                      (seq_samp[rec.id], rec.id)
                      ] = rec.seq

            samples.append(seq_samp[rec.id])

            if rec.id.split("___")[0] in seq_ids:
                seq_ids[seq_samp[rec.id]].append("%s___%s" %
                                                 (seq_samp[rec.id], rec.id))
            else:
                seq_ids[seq_samp[rec.id]] = ["%s___%s" %
                                             (seq_samp[rec.id], rec.id)]

        if len(set(samples)) == len(samples):
            with open("%s/%s.aln" % (outfolder, file_id), "w") as fout:
                for k in sequences:
                    fout.write(">%s\n%s\n" % (k, sequences[k]))
        else:
            sequences_df = pd.DataFrame.from_dict(sequences_df)
            clusters = paralogs_separator(seq_ids, sequences_df)
            for i, cluster in enumerate(clusters):
                with open("%s/%s_%d.aln" % (outfolder, file_id, i),
                          "w") as fout:
                    for k in cluster:
                        fout.write(">%s\n%s\n" % (k, sequences[k]))


if __name__ == '__main__':
    run()
