#!/usr/bin/env python

from os import system
from glob import glob
import click


@click.command()
@click.option(
    "-r",
    "--range",
    "crange",
    help="Range in for chr:start-end",
    required=True,
    default=None,
    type=str,
    show_default=True,
)
def main(crange):
    """
    Subrange a sorted BAM file
    """
    for fl in glob("*.bam"):

        system("samtools index {}".format(fl))
        system("samtools view -b -h -o tmp.bam {} {}".format(fl, crange))
        system("rm {}".format(fl))
        system("mv tmp.bam {}".format(fl))
        system("samtools index {}".format(fl))


if __name__ == "__main__":
    main()
