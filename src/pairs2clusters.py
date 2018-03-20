#!/usr/bin/env python
import click
import pandas as pd
import networkx as nx
from os import path

@click.command()
@click.option("-tbf", help="Tabfile containing two columns", type=str,
              default=None, show_default=True)
@click.option("-tbh", help="Table header", type=bool, default=True,
              show_default=True)
@click.option("-clf", help="Clusterfile Output", type=str, default=None,
              show_default=True)
def run(tbf, tbh, clf):
    if not tbf:
        exit("Input tab file is not given. Exiting . . . . .")
    if not clf:
        exit("Output cluster file is not given. Exiting . . . . .")
    if not path.isfile(tbf):
        exit("Given input file path either doesn't exist or it is not a file"
             " Exiting ....")
    # try:
    if tbh:
        data = pd.read_table(tbf)
    else:
        data = pd.read_table(tbf, header=None)
    graph = nx.Graph()
    graph.add_edges_from(data.values.tolist())
    connected_lists = nx.algorithms.components.connected.connected_components(graph)
#         for cn in connected_lists:
#             print(','.join(cn))
    with open(clf, "w") as fout:
        for lst in connected_lists:
            lt = ",".join(map(str, list(lst)))
            # print(lt)
            fout.write("%s\n" % lt)
#         # print("test")
    # except:
        # exit("Given file is not properly formatted. Exiting ....")

    # pass


if __name__ == '__main__':
    run()
