# TODO : Convert this file to see open frames in sequence using command line

import click
from Bio import SeqIO

@click.command()
@click.option("-infa", help="Input fasta", type=str, default=None,
              show_default = False)
@click.option("-ids", help="Display only selected sequences", type=str,
              default=None, show_default=True)
@click.option("-rc", help="Diplay reversed complement sequences", type=bool,
              default=False, show_default)
def run(infa, ids, rc):
    """Use less command in pipe mode."""
    assert infa, "Input fasta file not given. Exiting ..."
    try:
        # TODO: Put colours
        pass
    except:
        click.echo("Infile might not be in fasta or other option might not\
                   be valid. Exiting ...")

if __name__ == '__main__':
    run()
