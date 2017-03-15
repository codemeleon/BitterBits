import click
from Bio import SeqIO
from os import path


@click.command()
@click.option("-seqfile", help="Sequence file in genbank or fasta format",
              type=str, default=None, show_default=True)
@click.option("-seqfmt", help="Sequence format",
              type=click.Choice(['fasta', 'genbank']), default='fasta',
              show_default=True)
@click.option("-seqid", help="Sequence Id", type=str, default=None,
              show_default=True)
@click.option("-start", help="Start position (0 index based)",
              type=int, default=0, show_default=True)
@click.option("-end", help="End position (0 index based)",
              type=int, default=None, show_default=True)
@click.option("-rvcomp", help="Reverse complement (orly for DNA and RNA)",
              type=bool, default=False, show_default=True)
@click.option("--stdout", help="Print on screen", type=bool,
              default=True, show_default=True)
@click.option("-outfile", help="Output file. Valide when stdout==False",
              default=None, show_default=True)
def run(seqfile, seqfmt, seqid, start, end, rvcomp, stdout, outfile):
    """Extract genomic sequences for given range."""
    if not seqfile:
        click.echo("Sequence file not given. Exiting ....")
        exit(1)
    if not seqid:
        click.echo("Sequence id not given. Exiting ....")
        exit(1)
    if not start:
        click.echo("Start position not given. Considering from start of the\
                   sequences")
        start = 0
    if not end:
        click.echo("End position not given. Considering till end of the\
                   sequence")
    if end and end < start:
        click.echo("End has lower value than start. Exiting ....")
        exit(1)
    try:
        seq = SeqIO.to_dict(SeqIO.parse(seqfile, seqfmt))
        try:
            if not end:
                seq = seq[seqid].seq[start:]
            else:
                seq = seq[seqid].seq[start:end]
        except IOError:
            click.echo("Given sequence id is not in input file. Exiting")
            exit(1)
        if rvcomp:
            seq = seq.reverse_complement()
        if not stdout:
            with open(outfile, "w") as ofile:
                ofile.write(">%s\n%s\n" % (seqid, seq))
        else:
            click.echo(">%s\n%s\n" % (seqid, str(seq)))
    except:
        click.echo("Input file is not in given format. Exiting ....")
        exit(1)


if __name__ == '__main__':
    run()
