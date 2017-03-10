#!/usr/bin/env python
"""Select Uniq sequences."""

import click
from Bio import SeqIO
import pandas as pd


@click.command()
@click.option("-infile", help="Input fasta file", type=str, default=None,
              show_default=True)
@click.option("-outfile", help="Output fasta file", type=str, default=None,
              show_default=True)
def run(infile, outfile):
    """Generate file for unique sequences based on given input file."""
    assert infile, "Input file not fiven. Exiting ..."
    assert outfile, "Output file not fiven. Exiting ..."
    try:
        sequences = {"seq_id": [], "seq": []}
        for rec in SeqIO.parse(infile, "fasta"):
            sequences["seq_id"].append(rec.id)
            sequences["seq"].append(str(rec.seq))
        sequences = pd.DataFrame.from_dict(sequences)
        sequences = sequences.drop_duplicates(["seq"])
        with open(outfile, "w") as fout:
            for i, row in sequences.iterrows():
                fout.write(">%s\n%s\n" % (row["seq_id"], row["seq"]))
    except:
        click.echo("Enput files is not fasta or is empty")


if __name__ == '__main__':
    run()
