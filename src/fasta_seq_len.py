import click
from Bio import SeqIO

@click.command()
@click.option("--fasta_file", help="Input Fasta File", type=str,
              default=None, show_default=True)
def run(fasta_file):
    """Sequence length."""
    for rec in SeqIO.parse(fasta_file, "fasta"):
        print("%s\t%d"%(rec.id, len(rec.seq)))


if __name__ == '__main__':
    run()
