#!/usr/bin/env python

import click
from ete3 import Tree
from os import path

@click.command()
@click.option("-ti", help="Input newick tree", type=str, default=None,
              show_default=True)
@click.option("-to", help="Output newick Tree", type=str, default=None,
              show_default=True)
def run(ti, to):
    if not ti:
        exit("Input file not given. Exiting ...")
    if not to:
        exit("Output file not given. Exiting ....")
    if not path.isfile(ti):
        exit("Given path either doesn't exist or it is not a file. Exiting ..")
    try:
        tree = Tree(ti)
        tree.ladderize()
        tree.write(outfile=to)
    except:
        exit("Given file is not newick. Exiting . . .")


if __name__ == '__main__':
    run()

