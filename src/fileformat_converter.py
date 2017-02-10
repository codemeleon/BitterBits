import click
from Bio import SeqIO
from os import path


def gb2fasta(infile, outfile):
    """Genbank to Fasta."""
    with open(outfile, "w") as fout:
        for rec in SeqIO.parse(infile, "genbank"):
            fout.write(">%s\n%s\n" % (rec.id, rec.seq))
    return


@click.command()
@click.option("-infile", help="Input file", type=str,
              default=None, show_default=True)
@click.option("-infmt", help="Infile format", type=str,
              default='auto', show_default=True)
@click.option("-outfile", help="Output file path", type=str,
              default=None, show_default=True)
@click.option("-outfmt", help="Output file format", type=str,
              default='auto', show_default=True)
def run(infile, infmt, outfile, outfmt):
    """Convert file from one biological format to other."""
    if not infile:
        click.echo("Input file is not given. exiting ...")
        exit(1)
    if not path.isfile(infile):
        click.echo("Given file path dosn't exist. Exiting .....")
        exit(1)
    if not outfile:
        click.echo("Output file path not given")
        exit(1)
    if infmt:
        pass
    if outfmt:
        pass
    if infmt == 'genbank' and outfmt == 'fasta':
        gb2fasta(infile, outfile)


if __name__ == '__main__':
    run()
