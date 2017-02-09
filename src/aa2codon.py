import click
from os import path
from Bio import SeqIO


def amino2codon(aa, nuc):
    "Convert to Codon."
    codon = ''
    n = 0
    for a in aa:
        if a == '-':
            codon += '---'
        else:
            codon += nuc[3*n: 3*(n+1)]
            n += 1
    return codon


@click.command()
@click.option("--aaalnfile", help="Amino acid alignment file",
              type=str, default=None, show_default=True)
@click.option("--nucfile", help="Coding nucleotide folder",
              type=str, default=None, show_default=True)
@click.option("--codonalnfile", help="Codon alignment output",
              type=str, default=None, show_default=True)
def run(aaalnfile, nucfile, codonalnfile):
    """Code generate codon alignment base on amino acid alignments."""
    """ corresponding nucleotide sequences."""

    if aaalnfile is None or not path.isfile(aaalnfile):
        click.echo("Amino acid alignemt file is not given or path doesn't\
                   exist. Exiting ...")
        exit(1)
    elif nucfile is None or not path.isfile(nucfile):
        click.echo("Nucleotide folder is not given or path doesn't exist."
                   " Exiting ....")
        exit(1)
    elif not codonalnfile:
        click.echo("Codon file path not given. Exiting ...")
        exit(1)

    nuc_seq = {}
    for rec in SeqIO.parse(nucfile, 'fasta'):
        nuc_seq[rec.id] = rec.seq
    with open(codonalnfile, "w") as fout:
        for rec in SeqIO.parse(aaalnfile, "fasta"):
            fout.write(">%s\n%s\n" % (rec.id,
                                      amino2codon(rec.seq,
                                                  nuc_seq[rec.id])))


if __name__ == '__main__':
    run()
