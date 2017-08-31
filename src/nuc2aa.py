#!/usr/bin/env python

import click
from Bio import SeqIO
from Bio.Alphabet import IUPAC
import sys


@click.command()
@click.option("-inf", help="Nucleotide infile", type=str, default=None,
              show_default=True, required=True)
@click.option("-ouf", help="Amino acid infile. STDOUT if not given", type=str,
              default=None, show_default=True)
@click.option("-seqids", help="Comma separated seqids. All if not given",
              type=str, default=None, show_default=True)
@click.option("-frame", help="Translation frame",
              type=click.Choice(["1", "2", "3", "-1", "-2", "-3"]), default="1",
              show_default=True)
def run(inf, ouf, seqids, frame):
    """Converts Nulceotides file to Amino acid file."""
    if ouf:
        outfile = open(ouf, "w")
    else:
        outfile = sys.stdout
    frame = int(frame)
    rv = True if (frame< 0) else False
    frame = abs(frame) - 1
    if seqids:
        seqids = [s.strip() for s in seqids.split(',')]
    found = []
    for rec in SeqIO.parse(inf, 'fasta'):
        seq = rec.seq
        if rv:
            seq = rec.seq.reverse_complement()
        seq = seq[frame:]

        if seqids:
            if rec.id in seqids:
                found.append(rec.id)
                print(">%s\n%s\n" % (rec.id, seq.translate()), file=outfile)
        else:
            print(">%s\n%s\n" % (rec.id, seq.translate()), file=outfile)
    if seqids:
        notfound = ", ".join(set(seqids) -set(found))
        print(notfound, file=sys.stderr)

    if ouf:
        outfile.close()

if __name__ == '__main__':
    run()
