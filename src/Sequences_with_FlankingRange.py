import click
from Bio import SeqIO
from os import path, system
import pandas as pd

@click.command()
@click.option("-ref", help="Reference fasta genome", type=str,
              default=None, show_default=True)
@click.option("-query", help="Fasta Query File", type=str, default=None,
              show_default=True)
@click.option("-flank",  help="Fanking region length", type=int,
              default=0, show_default=True)
@click.option("-out", help="Output file", type=str, default=None,
              show_default=True)
@click.option("-blat", help="Blat path", type=str, default="blat",
              show_default=True)
def run(ref, query, flank, out, blat):
    """extract flanking flanking region squences from genome"""
    if not ref:
        click.echo("Ref file in not given. Exiting")
        exit(1)
    if not path.isfile(ref):
        click.echo("Given reference file doen't exits")
        exit(1)
    if not query:
        click.echo("Query file not given ....")
        exit(1)
    if not path.isfile(query):
        click.echo("Given query file path doesn't exist")
        exit(1)
    if flank < 0:
        click.echo("Negative Flank is not allowed. Exiting .....")
        exit(1)

    try:
        reference = SeqIO.to_dict(SeqIO.parse(ref, "fasta"))
    except TypeError:
        click.echo("Reference file is not fasta")
        exit(1)

    try:
        grbg = SeqIO.parse(query, "fasta")
    except TypeError:
        click.echo("Reference file is not fasta")
        exit(1)
    system("%s -noHead %s %s match.psl" % (blat, ref, query))
    data = pd.read_table("match.psl", header=None)
    data = data.sort_values([9, 0], ascending=[False, False])
    data = data.drop_duplicates([9])
    data = data[data[0] > 0.8*data[10]]
    outfile = open(out, "w")
    for i, row in data.iterrows():
        tseq = reference[row[13]].seq[row[15]-flank: row[16]+flank]
        if row[8] == '-':
            tseq = tseq.reverse_complement()
        outfile.write(">%s\n%s\n" % (row[9], tseq))
    outfile.close()
    system("rm match.psl")
    del grbg


if __name__ == '__main__':
    run()
