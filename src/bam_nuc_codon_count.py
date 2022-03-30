#!/usr/bin/env python

from genericpath import samefile
from glob import glob
from Bio import SeqIO, Seq
import sys
import pysam
import numpy as np
import pandas as pd
import click
from os import path

__author__ = "Anmol Kiran"
__organisation__ = (
    "Malawi-Liverpool-Wellcome Trust, Malawi; University of Liverpool, UK"
)
__github__ = "codemeleon"
__email__ = "akiran@mlw.mw"
__version__ = "0.0.1"


_codon_table = {
    "TAG": "*",
    "TAA": "*",
    "TGA": "*",
    "TTT": "F",
    "TTC": "F",
    "TTA": "L",
    "TTG": "L",
    "TCT": "S",
    "TCC": "S",
    "TCA": "S",
    "TCG": "S",
    "TAT": "Y",
    "TAC": "Y",
    "TGT": "C",
    "TGC": "C",
    "TGG": "W",
    "CTT": "L",
    "CTC": "L",
    "CTA": "L",
    "CTG": "L",
    "CCT": "P",
    "CCC": "P",
    "CCA": "P",
    "CCG": "P",
    "CAT": "H",
    "CAC": "H",
    "CAA": "Q",
    "CAG": "Q",
    "CGT": "R",
    "CGC": "R",
    "CGA": "R",
    "CGG": "R",
    "ATT": "I",
    "ATC": "I",
    "ATA": "I",
    "ATG": "M",
    "ACT": "T",
    "ACC": "T",
    "ACA": "T",
    "ACG": "T",
    "AAT": "N",
    "AAC": "N",
    "AAA": "K",
    "AAG": "K",
    "AGT": "S",
    "AGC": "S",
    "AGA": "R",
    "AGG": "R",
    "GTT": "V",
    "GTC": "V",
    "GTA": "V",
    "GTG": "V",
    "GCT": "A",
    "GCC": "A",
    "GCA": "A",
    "GCG": "A",
    "GAT": "D",
    "GAC": "D",
    "GAA": "E",
    "GAG": "E",
    "GGT": "G",
    "GGC": "G",
    "GGA": "G",
    "GGG": "G",
}


@click.command()
@click.option(
    "-bam",
    help="Bam files",
    default="/home/devil/Documents/Tools/BitterBits/src/test_data/"
    "K032623-rep-consensus_alignment_sorted.REF_NC_045512.2.bam",
    type=str,
    show_default=True,
)
@click.option(
    "-rid",
    help="Reference ID",
    type=str,
    default="NC_045512.2",
    show_default=True,
)
@click.option(
    "-ref",
    help="Reference fasta files",
    type=str,
    default="/home/devil/Documents/Tools/BitterBits/src/test_data/NC_045512.2.fasta",
    show_default=True,
)
@click.option(
    "-coor_range",
    help="Coordinates in the reference, zero index based, end exclusive",
    type=str,
    default="21000-25000",
    show_default=True,
)
@click.option(
    "-mut",
    help="Select locations with mismatches only",
    type=bool,
    default=False,
    show_default=True,
)
@click.option(
    "-gff",
    help="Gff Annotation File",
    type=str,
    default="/home/devil/Documents/Tools/BitterBits/src/test_data/genemap.gff",
    show_default=True,
)
@click.option(
    "--ignore_orphans",
    help="Ignore orphaned (Unpaired) reads",
    type=bool,
    default=False,
    show_default=True,
)
@click.option(
    "--min_mapping_quality",
    help="Mapping quality of reads",
    type=int,
    default=0,
    show_default=True,
)
@click.option(
    "--min_base_quality",
    help="Minimum base quality for correct base call",
    type=int,
    default=0,
    show_default=True,
)
@click.option(
    "--ignore_overlaps",
    help="Ignore paired overlapping reads",
    type=bool,
    default=False,
    show_default=True,
)
# TODO: Add samfiles relates conditions
@click.option("-codon", help="Select codon", type=bool, default=True, show_default=True)
def run(bam, rid, coor_range, mut, ref, gff, codon):
    """Expected to that bam file is sorted based on coordinate and indexed."""
    if not gff and codon:
        exit("Gff file is not given. Codon analysis is not possible.Exiting.")
    if gff:
        if not codon:
            print("No codon analysis. Ignoring gff file.")
        elif codon:
            if not path.exists(gff):
                exit("Gff file does not exist. Exiting.")
            elif not path.isfile(gff):
                print("Gff path is not a file")
            else:
                # NOTE: GFF to table
                with open(gff) as f:
                    gff_data = f.read()
                    gff_data = gff_data.split("##FASTA")[0]
                    gff_data = gff_data.split("\n")
                    gff_data = [x for x in gff_data if (
                        x != "" and x[0] != "#")]
                    gff_data = [x.split("\t") for x in gff_data]
                    gff_data = pd.DataFrame(
                        gff_data,
                        columns=[
                            "seq_id",
                            "source",
                            "feature",
                            "start",
                            "end",
                            "score",
                            "strand",
                            "frame",
                            "attribute",
                        ],
                    )
                    gff_data = gff_data.loc[gff_data["feature"] == "CDS"]
                    gff_data["start"] = gff_data["start"].astype(int)
                    gff_data["end"] = gff_data["end"].astype(int)
                    if rid not in gff_data["seq_id"].unique():
                        print("Reference sequence is not in gff file. Exiting.")
                        print("Reference", gff_data["seq_id"].unique())
                        exit()

            # TODO: Need details whether the fist and last are in the same CDS. If not, need to split the CDS.
            # TODO: If a region is not part of a CDS, consider that as part of intergenic region.
            print("Codon analysis is selected. Ignoring gff file.")

    sequences = {}
    for rec in SeqIO.parse(ref, "fasta"):
        sequences[rec.id] = rec.seq

    bam_files = []

    if not bam:
        print("Bam file not given.")
        sys.exit(0)
    if not path.exists(bam):
        print(f"Given bam file path {bam} doesn't exist.")
        sys.exit(0)
    if path.isdir(bam):
        bam_files = glob(f"{bam}/*.bam")
    if path.isfile(bam):
        bam_files = [bam]
        # print(f"Given bam path {bam} is not file")
        # sys.exit(0)
    print(bam)
    for bam in bam_files:
        # TODO: Check if bam file is sorted
        # TODO:check id bam file is indexed
        pass
    samfile = pysam.AlignmentFile(bam, "rb")
    # print(samfile.references)
    if rid not in samfile.references:
        print(f"Given reference {ref} not in given bam file")
        print("List of references")
        print(samfile.references)
        sys.exit(0)
    coor_range = coor_range.split("-")
    if len(coor_range) == 2:
        try:
            start = int(coor_range[0])
        except ValueError:
            exit("Range format is not correct, value before '-' is not numerical")
        try:
            end = int(coor_range[1])
        except Exception as e:
            exit("Range format is not correct, value after '-' is not numerical")
    elif len(coor_range) == 1:
        try:
            start = int(coor_range[0])
            end = int(coor_range[0])
        except:
            exit("Given coordinate is not numerical")
    else:
        print("Coordinate range is not in correct format")
        sys.exit(0)

    # TODO: Select the gene and it's regions
    # TODO: Check is coordinate is 0-based or 1-based
    # TODO: Define the codon position and amino acid position
    # TODO: Check if start and end are in the same gene. Might not be in the same due frame shift.
    # TODO: How will you handle the frame shift?
    # TODO: What if there are indels in the gene?
    # NOTE: Not for genome with splicing

    iter = samfile.pileup(
        rid,
        start,
        end + 1,
        ignore_orphans=False,
        min_mapping_quality=0,  # Play with all these three parameters
        min_base_quality=2,
        ignore_overlaps=False,
    )
    # TODO: First report the coordinates with mismatches only. Then report the codon and amino acid.
    # TODO: merge the nucletide which are part of the same codon.
    # TODO: Then check if they results in the same amino acid or not.
    mut = True
    coordinates_with_change = {}
    for pileupcol in iter:
        if (pileupcol.pos >= start) & (pileupcol.pos < end):
            # TODO: Include base call quality
            bases = {}
            for pread in pileupcol.pileups:
                if not pread.is_del and not pread.is_refskip:
                    if (
                        pread.alignment.query_sequence[pread.query_position]
                        not in bases
                    ):
                        bases[pread.alignment.query_sequence[pread.query_position]] = 0
                    bases[pread.alignment.query_sequence[pread.query_position]] += 1
            if mut:
                if len(bases) > 1:
                    coordinates_with_change[pileupcol.pos] = bases
                    if codon:
                        pass
                    else:
                        print(f"coordinate {pileupcol.pos}:", bases)
                else:
                    pass
            else:
                print(f"coordinate {pileupcol.pos}:", bases)

    # print(coordinates_with_change)
    final_table = {
        "Amino Acid Change": [],
        "Nucleotide Change": [],
        "Codon Change": [],
        "Samples": [],
        "Codon Count": [],
    }
    keys = set(coordinates_with_change)
    for _, row in gff_data.iterrows():
        selected_coordinates = keys & set(range(row["start"], row["end"]))
        if not selected_coordinates:
            continue
        else:
            # TODO: Add count of codon. Ignore with ambigious nucleotide and gaps
            # TODO: Check if there is any 3 base indels
            for selected_coordinate in selected_coordinates:
                shift = (selected_coordinate - row["start"] - 1) % 3
                iter = samfile.pileup(
                    rid,
                    selected_coordinate,
                    selected_coordinate + 1,
                    ignore_orphans=False,
                    min_mapping_quality=0,
                    min_base_quality=2,
                    ignore_overlaps=False,
                )
                codon_count = {}
                for pileupcol in iter:
                    if pileupcol.pos > selected_coordinate:
                        break
                    if pileupcol.pos < selected_coordinate:
                        continue
                    for pread in pileupcol.pileups:
                        if not pread.is_del and not pread.is_refskip:
                            codon = pread.alignment.query_sequence[
                                pread.query_position
                                - shift: pread.query_position
                                - shift
                                + 3
                            ]
                            if codon in _codon_table:
                                if codon not in codon_count:
                                    codon_count[codon] = 0
                                codon_count[codon] += 1

                            elif len(codon) > 3:
                                print(codon, selected_coordinate)
                            pass
                    # TODO: Perform reverse complement
                print(codon_count)
                row["strand"] = "-"
                if row["strand"] == "-":
                    codons = list(codon_count.keys())
                    # TODO: Generate conert and delete
                    for codon in codons:
                        codon_count[str(Seq.Seq(codon))] = codon_count[codon]
                    for codon in codons:
                        del codon_count[codon]
                    pass
                # TODO: Sort values in descreasing order
            # Coordinates to select
        pass
    final_table = pd.DataFrame(final_table)
    final_table = final_table.pivot_table(
        index=["Amino Acid Change", "Nucleotide Change", "Codon Change"],
        columns="Samples",
        values="Codon Count",
    )


if __name__ == "__main__":
    run()
