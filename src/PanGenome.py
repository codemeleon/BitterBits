#!/usr/bin/env python
"""Generate pangenome based on alignment files."""

import click
from Bio import AlignIO, SeqIO
from glob import glob
from os import path


@click.command()
@click.option("--alnfold", help="Fasta alignment folder", type=str,
              default=None, show_default=True)
@click.option("--outfile", help="Fasta pangenome alignment file",
              type=str, default=None, show_default=True)
@click.option("--samp_id_file", help="One id per line",
              type=str, default=None, show_default=True)
def run(alnfold, outfile, samp_id_file):
    r"""Generate pangenome based on given alignment files. Sequence ids are
    as samp___id."""
    assert alnfold, "Alignment folder not given. Exiting ...."
    assert outfile, "Output pangenome file not given. Exiting ...."
    assert path.isdir(alnfold), "Alignment folder doesn't exist. Exiting..."
    samples = []
    if samp_id_file and path.isfile(samp_id_file):
        with open(samp_id_file) as fin:
            for line in fin:
                samples.append(line[:-1])
        samples = set(samples)
    else:
        samples = set()
        for fl in glob("%s/*" % alnfold):
            seq_ids = SeqIO.to_dict(AlignIO.read(fl, "fasta")).keys()
            samp_ids = [seq_id.split("___")[0] for seq_id in seq_ids]
            samples |= set(samp_ids)
    genome_dict = {}
    for sample in samples:
        genome_dict[sample] = ""
    for fl in glob("%s/*" % alnfold):
        alignment = AlignIO.read(fl, "fasta")
        reported = []
        for rec in alignment:
            samp = rec.id.split('___')[0]
            reported.append(samp)
            genome_dict[samp] += str(rec.seq)

        for rec in samples - set(reported):
            genome_dict += '-' * alignment.get_alignment_length()

    with open(outfile, "w") as fout:
        for k in genome:
            fout.write(">%s\n%s\n" % (k, genome[k]))






if __name__ == '__main__':
    run()
