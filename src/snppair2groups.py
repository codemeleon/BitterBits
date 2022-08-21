#!/usr/bin/env python

# Requires input from pairsnp tool https://github.com/gtonkinhill/pairsnp

import pandas as pd
import click


@click.command()
@click.option('-inf',
              type=click.File('r'),
              required=True,
              help='Input file from pairsnp')
@click.option('-otf',
              type=click.File('w'),
              required=True,
              help='Output file to write to')
@click.option("-n",
              help="Number of SNPs to group together",
              default=10,
              type=int)
def run(inf, otf, n):
    """
    Convert snp pair file to groups file
    """
    # https://stackoverflow.com/questions/56567089/combining-lists-with-overlapping-elements
    df = pd.read_csv(inf, sep='\t')
    df.index = df.columns
    cells = [(df[col][df[col] <= 10].index[i], df.columns.get_loc(col))
             for col in df.columns
             for i in range(len(df[col][df[col] <= 10].index))]
    pooled = [set(subList) for subList in cells]
    merging = True
    while merging:
        merging = False
        for i, group in enumerate(pooled):
            merged = next((g for g in pooled[i + 1:] if g.intersection(group)),
                          None)
            if not merged: continue
            group.update(merged)
            pooled.remove(merged)
            merging = True
    otf.write("id\tgroup\n")
    for n, m in enumerate(pooled):
        for i in m:
            output_file.write(f'{i}\t{n}')


if __name__ == '__main__':
    run()
