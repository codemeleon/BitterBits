#!/usr/bin/env python

from Bio import SeqIO
import click
import pandas as pd
from os import system, makedirs


def gff2fna(gff, fna):
    """TODO: Docstring for gff2fna.

    :gff: TODO
    :fna: TODO
    :returns: TODO

    """
    pass


def gff2faa(gff, faa):
    """TODO: Docstring for gff2faa.
    :returns: TODO

    """
    pass


def faa2faauniq(faa, faa_u):
    """TODO: Docstring for faa2faauniq.
    :returns: TODO

    """
    pass


def faa2aln(faa, aln):
    """TODO: Docstring for faa2aln.

    :faa: TODO
    :aln: TODO
    :returns: TODO

    """
    pass


def faaaln2fnaaln(faa, fna):  # Need some additiona files
    """TODO: Docstring for faaaln2fnaaln.

    :faa: TODO
    :fna: TODO
    :returns: TODO

    """
    pass


@click.command()
@click.option(
    "-g",
    "--gff",
    "gff",
    help="Gff folder containing annotations and sequences per sample",
    type=str,
    default=None,
    show_default=True,
)
@click.option(
    "-a",
    "--aln",
    "aln",
    help="Alignment folder",
    type=str,
    default=None,
    show_default=True,
)
@click.option(
    "-c",
    "--cluster",
    "clst",
    help="Cluster file generated using roary",
    type=str,
    default=None,
    show_default=True,
)
def executer(gff, aln, clst):
    """TODO: Docstring for executer.

    :gff: folder containing gff files.
    :aln: alignent file for amino and nucleotides
    :clst: cluster file generated using roary
    :returns: None

    """
    pass


if __name__ == "__main__":
    executer()
