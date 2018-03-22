#!/usr/bin/env python

from Bio import Seq
import click
import pyperclip


def trans(seq, fm):
    if fm < 0:
        seq = seq.reverse_complement()
    seq = seq[abs(fm)-1:]
    print(">frame_%d\n%s\n"%(fm, seq.translate()))



@click.command()
@click.option("-fm", help="Frame (0 for all)",
              type=click.Choice([-3, -2, -1, 0, 1, 2, 3]),
              default=0, show_default=True)
def run(fm):
    """Translate clipboard nucleotide in different frame."""
    seq = Seq.Seq(pyperclip.paste().replace("\n","").replace("-",""))
    if not fm:
        for f in [1, 2, 3, -1, -2, 3]:
            trans(seq, f)
    else:
        trans(seq, fm)

if __name__=='__main__':
    run()
