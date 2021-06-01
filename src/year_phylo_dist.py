import click
from ete3 import Tree
import pandas as pd

@click.command()
@click.option("-t", help="Phylip tree file", type=str, default=None,
              show_default=True)
@click.option("-tab", help="Tab separated table containing sample in first \
              column and year(decimal) in the second", type=str, default=None,
              show_default=True)
@click.option("-ref", help="Reference in the Tree", type=str, default=None,
              show_default=True)
def run(t, tab, ref):
    "Calculates the distance of the of the samples with relative to given reference."
    table = pd.read_table(tab, header=None)
    leaf = list(table[0])
    leaf.append(ref)
    tree = Tree(t)
    t.prune(leaf)
    # TODO: Calculate the distance from the reference abd include in the table
    # TODO: Plot
    # Calculate the distance here


if __name__=='__main__':
    run()
