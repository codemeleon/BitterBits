import click
from ete3 import Tree
from os import path


@click.command()
@click.option(
    "-i",
    help="Input newick file having multiple leaves at same node",
    type=str,
    default=None,
    show_default=True,
)
@click.option(
    "-o",
    help="Ouput binary newick file",
    type=str,
    default=None,
    show_default=True,
)
def run(i, o):
    """This script generates binary tree to generate tangle gram using R
    package `dendextend`."""
    if not i:
        exit("Input file not given")
    if not path.exists(i):
        exit("Given input file path doesn't exist")
    if not path.isfile(i):
        exit("Given input path is not a file")
    if not o:
        exit("Output file path not given")
    try:
        tree = Tree(i)
        for n in tree.traverse("preorder"):
            if len(n.get_children()) > 2 and not n.is_root():
                leaf_names = n.get_leaf_names()
                for leaf in leaf_names[1:]:
                    dl = tree.search_nodes(name=leaf)[0]
                    dl.delete()
                tx = ""
                for x_ in leaf_names[1:-1]:
                    tx += f"({x_}:1e-9,"
                tx += (
                    f"{leaf_names[-1]}:1e-9"
                    + ")" * (len(leaf_names) - 2)
                    + ";"
                )
                tt = Tree(tx)
                n.add_child(tt)
        tree.write(outfile=o)
    except:
        print("Please check you input file format")


if __name__ == "__main__":
    run()
