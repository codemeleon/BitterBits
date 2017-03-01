#!/usr/bin/env python
"""Script to run multiAlign.py over a file."""

import click
from os import system, path
from multiprocessing import Pool
from functools import partial
from glob import glob


def runner(outfold, infile):
    """Worker."""
    file_base = path.split(infile)[1].split(".")[0]
    outfile = "%s/%s_aln.fasta" % (outfold, file_base)
    system("python multiAlign.py --i %s --o %s" % (infile, outfile))


@click.command()
@click.option("-inf", help="Input folder of fasta files",
              type=str, default=None, show_default=True)
@click.option("-outf", help="Output folder of fasta alignment",
              type=str, default=None, show_default=True)
@click.option("-ncor", help="Number of cores used in parallel",
              type=int, default=1, show_default=True)
def run(inf, outf, ncor):
    """Run multiAlign.py over a folder using multple core."""
    # TODO: add errors
    pool = Pool(ncor)
    func = partial(runner, outf)
    pool.map(func, glob("%s/*" % inf))


if __name__ == '__main__':
    run()
