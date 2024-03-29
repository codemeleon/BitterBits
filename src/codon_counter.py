#!/usr/bin/env python

import click
from os import path, system
import pysam as ps
import pandas as pd
import numpy as np
from Bio import SeqIO


def check_reference(mut_file, ref_file, gff_file, bam_file):
    """Check if reference name is same in all the input files."""
    fasta_ref = []
    gff_ref = []
    bam_ref = []
    for rec in SeqIO.parse(ref_file, "fasta"):
        fasta_ref.append(rec.id)

    fasta_ref = set(fasta_ref)

    gff = pd.read_table(gff_file, comment="#")
    gff_ref = set(gff[0])

    mut_ref = []
    with open(mut_file) as fin:
        for line in fin:
            ref = line.split(":")[0]
            mut_ref.append(ref)
    mut_ref = set(mut_ref)
    bamfile = ps.AlignmentFile(bam_file, "rb")
    bam_ref = set(bamfile.references)
    to_exit = 0

    additiona_mut_ref = mut_ref - fasta_ref
    if len(additiona_mut_ref):
        print(additiona_mut_ref, "are not present in given reference fasta file")
        to_exit = 1

    additiona_mut_ref = mut_ref - gff_ref
    if len(additiona_mut_ref):
        print(additiona_mut_ref, "are not present in given gff file")
        to_exit = 1

    additiona_mut_ref = mut_ref - bam_ref
    if len(additiona_mut_ref):
        print(additiona_mut_ref, "are not present in given bam file")
        to_exit = 1

    return to_exit


def name_the_gene(mut_file, gff_file):
    """TODO: Docstring for name_the_gene.

    :mut_file: TODO
    :gff_file: TODO
    :returns: TODO

    """
    mut_loc = {"ref": [], "pos": []}
    with open(mut_file) as fin:
        for line in fin:
            ref, pos = line[:-1].split(":")
            pos = int(pos)
            mut_loc["ref"].append(ref)
            mut_loc["pos"].append(pos)
    mut_loc = pd.DataFrame(mut_loc)
    gff = pd.read_table(gff_file, comment="#")
    gff["gene"] = ""  # TODO: Report the gene name in the table
    gene_set = []
    for _, row in mut_loc.iterrows():
        genes = set(
            gff.loc[
                (gff["ref"] == row["ref"])
                & (gff["start"] <= row["pos"])
                & (gff["end"] > row["pos"]),
                "gene",
            ].values
        )
        gene_set.append(genes)
    mut_loc["gene"] = gene_set
    # TODO instead of list of genes, convert them to list of coordinates

    return mut_loc, gff


def check_cds(mut_file, gff_file):
    """TODO: Docstring for check_cds.

    :mut_file: TODO
    :gff_file: TODO
    :returns: TODO

    """
    gff = pd.read_table(gff_file, comment="#")
    gff = gff[gff[10] == "CDS"]

    not_cds = 0
    with open(mut_file) as fin:
        for line in fin:
            chrom, coor = line[:-1].split(":")
            coor = int(coor)
            t_gff = gff[(gff[0] == chrom) & (
                gff[10] <= coor) & (gff[11] > coor)]
            if not len(t_gff):
                print(f"{line[:-1]} is not part of CDS.")
                not_cds = 1
    return not_cds


def codon_change_distributions(mut_file, gff_file, ref_file, bam_file):
    """TODO: Docstring for codon_change_distributions.

    :mut_file: TODO
    :gff_file: TODO
    :ref_file: TODO
    :bam_file: TODO
    :returns: TODO

    """
    bamfile = ps.AlignmentFile(bam_file, "rb")
    mut_loc, gff = name_the_gene(mut_file, gff_file)
    ref_seq = {}
    for rec in SeqIO.parse(ref_file, "fasta"):
        ref_seq[rec.id] = rec.seq

    for _, row in mut_loc.iterrows():
        start, end, strand = gff.head()  # TODO:Change it to relevent values
        gene_seq = ref_seq[row["gene"]][start:end]
        if strand == "-":
            gene_seq = gene_seq.reverse_complement()

        gene_seq_prot = gene_seq.translate()
        # CAUTION: beware of nucleotide the the begining of the reads in case of gene is on reverse strand
        # TODO: Provide ref id, gene id, position in ref, position in gene, stand of the gene on reference

        for gene in row["genes"]:

            pass
        pass

    # TODO: Start with gff file
    # TODO: Indenfify coorindtaes which can be groupped together
    # TODO: Report the changes
    with open(mut_file) as fin:
        for line in fin:
            chrom, coor = line[:-1].split(":")
            coor = int(coor)
            # TODO: check whether nucleotide are part of same codon
            # NOTE: Considering no frameshift in coding region and is unique

            pass

    pass


@click.command()
@click.option(
    "--mut",
    "-m",
    "mut_file",
    help="Mutation position files. 0 index based",
    type=str,
    default=None,
    show_default=True,
)
@click.option(
    "--ref",
    "-r",
    "ref_file",
    help="Reference fasta file.",
    type=str,
    default=None,
    show_default=True,
)
@click.option(
    "--gff",
    "-g",
    "gff_file",
    help="Gff file consisting CDS info",
    type=str,
    default=None,
    show_default=True,
)
@click.option(
    "--bam",
    "-b",
    "bam_file",
    help="Bam file, generated by mapping short reads against a reference",
    type=str,
    default=None,
    show_default=True,
)
def run(mut_file, ref_file, gff_file, bam_file):
    """TODO: Docstring for ru.

    :mut_file: TODO
    :returns: TODO

    """
    if check_reference(mut_file, ref_file, gff_file, bam_file):
        exit("Exiting. . . .")
    if check_cds(mut_file, gff_file):
        exit("Exiting . . . .")

    if not path.isfile(f"{ref_file}.fai"):
        system(f"samtools faidx {ref_file}")

    # TODO: Check whether bam file is sorted.
    # TODO: Sort given coordinates in mut file.
    # TODO: Check whether multiple coordinates are part of same codon

    pass


if __name__ == "__main__":
    run()
