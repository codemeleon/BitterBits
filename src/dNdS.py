#!/usr/bin/env python


import click
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from itertools import product
import networkx as nx
from fractions import Fraction
from functools import partial
from multiprocessing import Pool

def codon_pair_dnds_dict():
    # based on
    # http://bioinformatics.cvr.ac.uk/blog/calculating-dnds-for-ngs-datasets/
    graph = nx.Graph()
    codons = list(map(''.join, product("ATGC", repeat=3)))
    for c1 in codons:
        for c2 in codons:
            if c1 <= c2:
                continue
            if np.sum(np.array(list(c1)) != np.array(list(c2))) == 1:
                graph.add_edge(c1, c2)
    nd = {}
    sd = {}
    for c1 in codons:
        if c1 not in nd:
            nd[c1] = {}
            sd[c1] = {}
        for c2 in codons:
            if c2 not in nd[c1]:
                nd[c1][c2] = 0
                sd[c1][c2] = 0
            shortest_paths = list(nx.all_shortest_paths(graph, source=c1,
                                                        target=c2))
            dn = 0
            ds = 0
            for shortest_path in shortest_paths:
                # print(shortest_path)
                for pos in range(len(shortest_path)-1):
                    if (Seq(shortest_path[pos]).translate() ==
                            Seq(shortest_path[pos + 1]).translate()):
                        ds += 1
                    else:
                        dn += 1
            # if c1=="ATG" and c2=="TAG":
            #     print(dn, ds, len(shortest_paths))
            nd[c1][c2] = Fraction(dn, len(shortest_paths))
            sd[c1][c2] = Fraction(ds, len(shortest_paths))
    S = {}
    N = {}
    for c1 in codons:
        if c1 not in S:
            S[c1] = 0
        for c2 in graph[c1]:
            if Seq(c1).translate() == Seq(c2).translate():
                S[c1] += Fraction(1, 3)
        N[c1] = 3 - S[c1]
    return {'nd': nd, 'sd': sd, 'S': S, 'N': N, 'graph':graph}

# x = codon_pair_dnds_dict()
# x['nd']["ATG"]["TAG"]
# x['graph']["ATG"]["TAG"]
# list(nx.all_shortest_paths(x['graph'], source='ATG', target='ATG'))

def calculate(cpdd, seq1, seq2):
    # Senquence must be of same length and length must be mutiple of 3
    if len(seq1) != len(seq2) or len(seq1) % 3 != 0:
        exit(0)
    seq1 = [seq1[i:i+3] for i in range(0, len(seq1), 3)]
    seq2 = [seq2[i:i+3] for i in range(0, len(seq2), 3)]
    # print(seq2)
    # print(len(seq1), len(seq2))
    pair = pd.DataFrame.from_dict({'seq1': seq1, 'seq2': seq2})
    pair = pair[~((pair["seq1"] == '---') | (pair["seq2"] == '---'))]
    pair = pair[pair["seq1"] != pair["seq2"]]
    nd, sd, S, N = 0, 0, 0, 0
    if len(pair):
        for i, row in pair.iterrows():
            if ('-' in row["seq1"]) or ('-' in row["seq2"]):
                print("Check you alignment.\
                      There is a frame shift in alignment.")
                exit(1)
            # print(row["seq1"], row["seq2"])
            # print(cpdd['nd'][row['seq1']][row['seq2']], cpdd['sd'][row['seq1']][row['seq2']], "Anmol")
            S += cpdd['S'][row['seq1']]
            N += cpdd['N'][row['seq1']]
            nd += cpdd['nd'][row['seq1']][row['seq2']]
            sd += cpdd['sd'][row['seq1']][row['seq2']]
        # print(nd*1.0/ N, sd*1.0/ S)
            # print(cpdd['sd'][row['seq1']][row['seq2']])
        # print(float(1 - (sd / S) * Fraction(4, 3)))
            # print(sd,  S)
        S = float(S)
        N = float(N)
        nd = float(nd)
        sd = float(sd)
    return sd, nd, S, N
        # if (sd/S) == 0 and (nd/N) == 0:
        #     return 1
        # if (sd/S) == 0:
        #     return 1e6
        # dnds = (nd/N)/(sd/S)
        # return dnds
        # print(S, N, nd, sd, (nd*4) / (3*N), (sd*4) / (3*S))
        # pN = np.log(1 - (nd*4) / (3*N))
        # pS = np.log(1 - (sd*4) / (3*S))
    #     # print(pN, pS)
    #     if pS == 0 and pN == 0:
    #         # print("anmol")
    #         dnds = 1
    #     elif pS == 0:
    #         dnds = 1e6
    #     elif pN == 0:
    #         dnds = 0
    #     else:
    #         dnds = float(pN / pS)
    # else:
    #     dnds = 1
    # return 1


@click.command()
@click.option("-inf", help="Input fasta file", type=str, default="/home/devil/Documents/Jen/GPS/Fol/DNDS/tmp.fa",
              show_default=True)
@click.option("-cor", help="Number of CPU to use", type=int, default=1,
              show_default=True)
@click.option("-plot", help="Heatmap plot", type=bool, default=False,
              show_default=True)
def run(inf, cor, plot):
    if not inf:
        click.echo("Input file not given. Exiting")
        exit(1)
    try:
        sequencesgen = SeqIO.parse(inf, "fasta")
    except IOError:
        click.echo("Given file is not fasta. Exiting.")
        exit(1)
    sequences = {}
    for k in sequencesgen:
        sequences[k.id] = str(k.seq).upper()
    index = list(sequences.keys())
    cpdd = codon_pair_dnds_dict()
    pool = Pool(cor)
    dnds_dict = {}
    # print(list(sequences.values()))
    for k1 in sequences:
        func = partial(calculate, cpdd, sequences[k1])
        # print(pool.map(func, list(sequences.values())))
        dnds_dict[k1] = list(pool.map(func, list(sequences.values())))
    dnds_dict = pd.DataFrame.from_dict(dnds_dict)
    dnds_dict.index = index
    dnds_dict.to_csv("dnds_tab.tsv", sep="\t")


if __name__ == '__main__':
    run()
