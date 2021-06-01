#!/usr/bin/env python



import click
import pandas as pd
import numpy as np
from os import system, path




def transform(i, o):
    """Generate organised file based on the largest sequences

    :i: TODO
    :o: TODO
    :returns: TODO

    """
    system(f"blat {i} {i} {i}.psl")
    data = pd.read_table(f"{i}.psl", header=None)
    system(f"rm {i}.psl")
    data = data[data[9]!=data[13]]
    data = data.sort_values(0,
            ascending=False)
    data = data.drop_duplicates([9,13])
    data.loc[data[9]<data[13],[9,13] ] =  data.loc[data[9]<data[13],
            [13, 9] ].values
    data = data.drop_duplicates([9,13])
    # TODO: Get the longest sequence
    # TODO: Store the sequences in dictionary
    data = data[data[9]==ids[0]]
    with open(o, "w") as fout:
        fout.write(f">{}\n{}\n")
        for _, row in data.iterrows():
            pass














@click.command()
@click.option("-i",
      help="Input fasta",
      type=str,
      default=None,
      show_default=True)
@click.option(
    "-o",
    help="Output file",
    type=str,
    default=None,
    show_default=True)
def run(i, o):
    """This script organise the sequence based on the largest sequence."""
    pass

if __name__=='__main__':
    run()
