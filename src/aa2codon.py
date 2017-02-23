#! /usr/bin/env python
"""Protein to Codon."""
import click
from glob import glob
from os import path, makedirs
from Bio import SeqIO


def amino2codon(aa, nuc):
    """Convert to Codon."""
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
@click.option("--aaalnfold", help="Amino acid alignment folder",
              type=str, default=None, show_default=True)
@click.option("--nucfold", help="Coding nucleotide folder",
              type=str, default=None, show_default=True)
@click.option("--cdsalnfold", help="Codon alignment output folder",
              type=str, default=None, show_default=True)
@click.option("--ofilext", default="Output File extension",
              type=str, default="cds", show_default=True)
def run(aaalnfold, nucfold, cdsalnfold, ofilext):
    """Code generate codon alignment base on amino acid alignments."""
    """ corresponding nucleotide sequences."""

    assert aaalnfold, "Animo acid aligment folder not give. Exiting .."
    assert nucfold, "DNA input folder not give. Exiting .."
    assert cdsalnfold, "Codon Alignment output folder not give. Exiting .."
    assert path.isdir(aaalnfold), "Given amino acid alignment folder doesn't\
        exist. Exiting ..."
    assert path.isdir(nucfold), "Given nucleotide folder doesn't\
        exist. Exiting ..."
    assert not makedirs(cdsalnfold, exist_ok=True), "Unable to create output\
        folder. Please check permission. Exiting ..."

    nuc_seq = {}
    for fl in glob("%s/*" % nucfold):
        for rec in SeqIO.parse(fl, 'fasta'):
            nuc_seq[rec.id] = rec.seq

    for fl in glob("%s/*" % aaalnfold):
        file_id = path.split(fl)[1].split(".")[0]
        with open("%s/%s.%s" % (cdsalnfold, file_id, ofilext), "w") as fout:
            for rec in SeqIO.parse(fl, "fasta"):
                fout.write(">%s\n%s\n" %
                           (rec.id,
                            amino2codon(rec.seq,
                                        nuc_seq[rec.id])))


if __name__ == '__main__':
    run()
