#!/usr/bin/env python

import click
from Bio import SeqIO, Seq

def cds_aa(seq):
    frame, mx = 0,0
    for i in range(3):
        translated = str(seq[i:].translate())
        translated = translated.split("*")
        mx_len = max(map(len, translated))
        if mx_len > mx:
            mx = mx_len
            frame = i
    seq = seq[frame:]
    translated = str(seq.translate()).split("*")
    len_list = list(map(len, translated))
    mx_len = max(len_list)
    for i,lgnth in enumerate(len_list):
        if lgnth == mx_len:
            break
    # TODO: Extract the cds
    # TODO: Translate it


    return cds, aa


@click.command()
@click.option("-infa", help="Input nucleotide fasta", type=str, default=None,
              show_default=True)
@click.option("-ocds", help="Output CDS fasta", type=str, default=None,
              show_default=True)
@click.option("-ofaa", help="Output Amino Acids Fasta", type=str, default=None,
              show_default=True)
def run(infa, ofna, ofaa):
    """Extract longest Amino acid sequence and corresponing codon for
    Positive strand. It can create trouble is two open reading frames are
    of the same length. Use it with caution"""
    try:
        sequence = SeqIO.parse(infa,"fasta")
    except:
        exit("Either input nucleotide file not given or file does exists or \
             given file is not fasta. Exiting . . . . .")
    seq_out = {}
    for rec in sequence:
        seq_out[rec.id] = Seq.Seq(str(rec.seq).replace("-","").upper())
    del sequence
    seq_cds = {}
    seq_aa = {}
    for k in seq_out:




if __name__='__main__':
    run()
