#!/usr/binenv python

from Bio import SeqIO
import click
from os import path


@click.command()
@click.option("--inputgb", help="Input genbank file", type=str, default=None,
              show_default=True)
@click.option("--outputfasta", help="Output fasta file", type=str,
              default=None, show_default=True)
@click.option("--seq", help="Type of sequence",
              type=click.Choice(["genome","gene","CDS","protein"]),
              default='genome', show_default=True)
def run(inputdb, outputfasta, seq):
    "Extract sequence from genbank file"
    if not inputdb:
        click.echo("Input genbank file not given. Exiting ....")
        exit(1)
    elif not outputfasta:
        click.echo("Output fasta file not given. Exiting ....")
        exit(1)
    else:
        pass
    try:
        if seq == 'genome':
            SeqIO.convert(inputgb, "genbank", outputfasta, "fasta")
        else:
            record = SeqIO.read(inputgb, "genbank")
            for feature in record.features:
                if feature.type in ["gene", "CDS"]:
                    seq = record.seq[feature.location.start:
                                     feature.location.end]
                    if feature.location.strand == -1:
                        seq = seq.reverse_complement()
                    seq = str(seq).upper()
                    if feature.type == "CDS" and seq == "protein":
                        outputfasta.write(">%s\t%s\t%s\n%s\n\n" %
                                          (feature.qualifiers["locus_tag"][0],
                                           feature.qualifiers["gene"][0] if
                                           "gene" in feature.qualifiers.keys()
                                           else "",
                                           feature.qualifiers["product"][0] if
                                           "product" in
                                           feature.qualifiers.keys()
                                           else "",
                                           str(feature.qualifiers[
                                               "translation"][0])))
                    elif feature.type == seq:
                        outputfasta.write(">%s\t%s\t%s\n%s\n\n" %
                                          (feature.qualifiers["locus_tag"][0],
                                           feature.qualifiers["gene"][0] if
                                           "gene" in feature.qualifiers.keys()
                                           else "",
                                           feature.qualifiers["product"][0] if
                                           "product" in
                                           feature.qualifiers.keys() else "",
                                           seq))
    except:
        pass
        # Do something silly


if __name__ == '__main__':
    run()
