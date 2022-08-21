import click
import pandas as pd
from Bio import SeqIO
import swifter


def most_common(lst):
    current_count = 0
    current_base = ''
    for x in 'ATGC':
        count = list(lst).count(x)
        if count > current_count:
            current_count = count
            current_base = x
    return current_base


@click.command()
@click.option(
    "-i",
    "--in",
    "inputfile",
    help="Fasta input file",
    type=str,
    default=None,
    show_default=True,
)
@click.option(
    "-o",
    "--out",
    "outputfile",
    help="Fasta output file",
    type=str,
    default=None,
    show_default=True,
)
def run(inputfile, outputfile):
    """Generates concensus sequence

    :inputfile: Input fasta file
    :outputfile: Output fasta file
    :returns: None
    """
    sequences = {}
    # n = 0
    for rec in SeqIO.parse(inputfile, "fasta"):
        sequences[rec.id] = list(str(rec.seq).upper())
        # n += 1
        # if n == 10:
        # break
    sequences = pd.DataFrame(sequences)
    sequences = ''.join(sequences.swifter.apply(most_common, axis=1).values)
    with open(outputfile, 'w') as fout:
        fout.write(f">seq\n{sequences}\n")


if __name__ == "__main__":
    run()
