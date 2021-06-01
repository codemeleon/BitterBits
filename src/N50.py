from glob import glob
import click
from os import path
from Bio import SeqIO
import numpy as np


@click.command()
@click.option("-f", help="Fasta file of containing folder", type=str, default=None,show_default=True)
@click.option("-o", help="Output csv file", type=str, default=None,show_default=True)
def run(f, o):
    """Generates N50 from Fasta files"""
    if not f:
        exit("Input file/folder path not given. Exiting. . . . ")
    elif not o:
        exit("Output file path not given. Exiting. . . . ")
    elif not path.exists(f):
        exit("Given input path does not exist. Exiting. . . . ")
    if path.isfile(f):
        try:
            seq_len = []
            for rec in SeqIO.parse(f, "fasta"):
                seq_len.append(len(rec.seq))
            seq_len = np.cumsum(np.sort(seq_len)[::-1])
            N50 = np.sum(seq_len<seq_len[-1]/2.)+1
            with open(o, "w") as fout:
                fout.write(f"{f}\t{N50}\t{seq_len[-1]}\n")
        except:
            print(f"{f} is not Fasta file")
    else:
        with open(o, "w") as fout:
            for fl in glob(f"{f}/*"):
                try:
                    seq_len = []
                    for rec in SeqIO.parse(fl, "fasta"):
                        seq_len.append(len(rec.seq))
                    seq_len = np.cumsum(np.sort(seq_len)[::-1])
                    N50 = np.sum(seq_len<seq_len[-1]/2.)+1
                    fout.write(f"{fl}\t{N50}\t{seq_len[-1]}\n")
                except:
                    print(f"{fl} is not Fasta file")

if __name__=='__main__':
    run()