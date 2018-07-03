#!/usr/bin/env python
import click
import pandas as pd
from Bio import SeqIO
from os import path
from ete3 import Tree

@click.command()
@click.option("-intree", help="Input Newick Tree", type=str, default=None,
              show_default=True)
@click.option("-aln", help="Corresponing alignment fasta file", type=str,
              default=None, show_default=True)
@click.option("-outtree", help="Output Tree", type=str, default=None,
              show_default=True)
def run(intree, aln, outtree):
    """This script includes additional leaves nodes which have been removed by
    RAxML due to identical sequences"""
    if not intree:
        exit("Input tree not given")
    elif not path.isfile(intree):
        exit("Given input tree is not a file or path doesn't exists")
    if not aln:
        exit("Input alignment not given")
    elif not path.isfile(aln):
        exit("Given input alignment is not a file or path doesn't exists")
    if not outtree:
        exit("Output tree not given")

    try:
        tree = Tree(intree)
    except:
        exit("Gven file is not phylip")
    try:
        sequences = SeqIO.parse(aln, "fasta")
    except:
        exit("Give file is not fasta")

    sequences_df = {'Id':[], "seq":[]}
    for k in sequences:
        sequences_df['Id'].append(k.id)
        sequences_df['seq'].append(str(k.seq).upper())
    sequences_df = pd.DataFrame(sequences_df)
    sequences_groups = sequences_df.groupby("seq")["Id"].apply(list)\
        .reset_index()["Id"]
    sequences_groups = sequences_groups[sequences_groups.apply(len)>1]
    sequences_groups = list(sequences_groups)
    for leaf in tree.get_leaves():
        for lst in sequences_groups:
            if leaf.name in lst:
                for l in lst:
                    if l == leaf.name:
                        continue
                    leaf.add_sister(name=l,dist=0)
                sequences_groups.remove(lst)
    tree.write(outfile=outtree)



if __name__=='__main__':
    run()
