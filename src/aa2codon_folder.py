import click
from glob import glob
from os import path, makedirs
from multiprocessing import Pool
import pandas as pd
from Bio import SeqIO

def amino2codon(aa, nuc):
    codon = ''
    n = 0
    for a in aa:
        if a == '-':
            codon += '---'
        else:
            codon += nuc[3*n: 3(n+1)]
            n += 1
    return codon

def nuc_table(nucfold):
    nuc_dict = {}
    for fl in glob("%s/*" % nucfold):
        for rec in SeqIO.parse(fl, 'fasta'):
            nuc_dict[rec.id] = rec.seq

@click.command()
@click.option("--aaalnfold", help="Amino acid alignment folder'alnout'",
              type=str, default=None, show_default=True)
@click.option("--nucfold", help="Coding nucleotide folder",
              type=str, default=None, show_default=True)
@click.option("--codonalnfold", help="Codon alignment output\
              [default: aaanfold]", type=str, default=None)
def run(aaalnfold, nucfold, codonalnfold):
    """Code generate codon alignment base on amino acid alignments and"""
    """ corresponding nucleotide sequences."""

    if aaalnfold is None or not path.isdir(aaalnfold):
        click.echo("Amino acid alignemt folder is not given or path doesn't\
                   exist")
        exit(1)
    elif nucfold is None or not path.isdir(nucfold):
        click.echo("Nucleotide folder is not given or path doesn't exist")
        exit(1)
    elif not codonfolder:
        codonfolder = aaalnfold
    elif not path.isdir(codonalnfold):
        click.echo("Codon folder doesn't exist. Creating")
        makedirs(codonalnfold)

    nuc_seq = nuc_table(nucfold)
    for fl in glob("%s/*" % aaalnfold):
        outfile = open("%s/%s.codon" %
                       (codonalnfold,
                        path.split(fl)[1].split(".fa")[0]), "w")
        for rec in SeqIO.parse(fl, "fasta"):
            outfile.write(">%s\n%s\n" % (rec.id,
                                         amino2codon(rec.seq,
                                                     nuc_seq[rec.id])))
        outfile.close()


if __name__ == '__main__':
    run()
