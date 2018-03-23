#!/usr/bin/env python

import click
import pyperclip
from os import path
import pandas as pd
from Bio import AlignIO

@click.command()
@click.option("-inf", help="Alignment Fasta", type=str, default=None,
               show_default=True)
@click.option("-otf", help="Output SNP alignment Fasta", type=str,
               default=None, show_default=True)
def run(inf, otf):
    """Outputs alignment with variable location from the given alignment

    If file name is in clipboard:

        if out file not given:

            if want to print on screen:

                python this_script

            else:

                python this_script > outputfile

        else:

            python -otf outputfile

    If input file name:

        if out file not given:

            if want to print on screen:

                python this_script -inf inputfile

            else:

                python this_script -inf inputfile > outputfile

        else:

            python -inf inputfile -otf outputfile


    """
    if not inf:
        inf = pyperclip.paste()
    if not path.isfile(inf):
        exit("Given path or clipbord is not valid file path. Exiting.....")
    try:
        alignment = {}
        for rec in AlignIO.parse(inf, "fasta"):
            alignment[rec.id] = list(rec.seq)
        alignment = pd.DataFram(alignment)
        alignment = alignment[alignment.apply(lambda x: len(set(x)))!=1]
        if not otf:
            with open(otf, "w") as fout:
                for k in alignment.columns:
                    fout.write(">%s\n%s\n" % (k, ''.join(alignment[k])))
            pass
        else:
            for k in alignment.columns:
                print(">%s\n%s\n" % (k, ''.join(alignment[k])))
    except:
        exit("Probably given file is not fasta or not of sequencs in "
             "alignment of not same length. Exiting......")



if __name__=='__main__':
    run()
