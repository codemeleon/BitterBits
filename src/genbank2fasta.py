##!/usr/binenv python

from Bio import SeqIO
import click
from os import path


@click.command()
@click.option("--inputgb", help="Input genbank file", type=str, default=None,
              show_default=True)
@click.option("--outputfasta", help="Output fasta file", type=str,
              default=None, show_default=True)
@click.option("--seqtype", help="Type of sequence",
              type=click.Choice(["genome", "gene", "CDS", "protein"]),
              default='genome', show_default=True)
@click.option("--seqid", help="Sequence id", type=click.Choice(["protein_id"]),
              default="protein_id", show_default=True)
@click.option("--fileprefix", help="Add file name separated iwth ___ in\
              sequence name", type=bool, default=False, show_default=True)
def run(inputgb, outputfasta, seqtype, seqid, fileprefix):
    "Extract sequence from genbank file"
    assert inputgb, "Input genbank file not given. Exiting ...."
    assert outputfasta, "Output fasta file not given. Exiting ...."

    # try:
    if seqtype == 'genome':
        SeqIO.convert(inputgb, "genbank", outputfasta, "fasta")
    else:
        outputfasta = open(outputfasta, "w")
        input_file_id = path.split(inputgb)[1].split('.')[0]
        records = SeqIO.parse(inputgb, "genbank")
        for record in records:
            for feature in record.features:
                # print("Anmol", feature.type)
                if feature.type in ["gene", "CDS"]:
                    seq = record.seq[feature.location.start:
                                     feature.location.end]
                    if feature.location.strand == -1:
                        seq = seq.reverse_complement()
                    seq = str(seq).upper()
                    # if feature.type == "CDS":
                    if ((seqtype in ['CDS', 'protein']) and
                        (('translation' not in feature.qualifiers.keys()) or
                         ('protein_id' not in feature.qualifiers.keys()))):
                        continue

                    if feature.type == "CDS":
                        if seqtype == "protein":
                            if fileprefix:
                                outputfasta.write(">%s___%s\n%s\n\n" % (
                                    input_file_id,
                                    feature.qualifiers["protein_id"][0],
                                    str(feature.qualifiers["translation"][0])))
                            else:
                                outputfasta.write(">%s\n%s\n\n" % (
                                    feature.qualifiers["protein_id"][0],
                                    str(feature.qualifiers["translation"][0])))

                        elif seqtype == "CDS":
                            if fileprefix:
                                outputfasta.write(">%s___%s\n%s\n\n" % (
                                    input_file_id,
                                    feature.qualifiers["protein_id"][0],
                                    seq))
                            else:
                                outputfasta.write(">%s\n%s\n\n" % (
                                    feature.qualifiers["protein_id"][0],
                                    seq))
        outputfasta.close()


if __name__ == '__main__':
    run()
