#!/usr/bin/env python
"""Generate table for genes in blast search (table output)."""
import click
from glob import glob
import pandas as pd
from os import path


@click.command()
@click.option("-inf", help="Input folder", type=str, default=None,
              show_default=True)
@click.option("-ext", help="File extention to select files", type=str,
              default="blast8", show_default=True)
@click.option("-idtt", help="Identity relative to query", type=float,
              default=0.5, show_default=True)
@click.option("-outfile", help="Output file path", type=str,
              default=None, show_default=True)
def run(inf, ext, idtt, outfile):
    """Generate table for genes in blast search."""
    # TODO: Code need some improvements. Make it work later
    assert inf, "Input folder not given. Exiting ...."
    assert path.isdir(inf), "Given folder path is not a folder. Exiting ..."
    assert outfile, "Output file not given. Exiting ...."
    final_table = pd.DataFrame()
    for fl in glob("%s/*.%s" % (inf, ext)):
        table = pd.read_table(fl, header=None)
        table = table[table[2] > (idtt * 100.)][[0, 1]]
        if not len(table):
            continue
        table = table.groupby(1)[0].apply(','.join).reset_index()
        # print(table)
        table.index = table[1]
        del table[1]
        table = table.T
        table["sample"] = [path.split(fl)[1].split(".%s" % ext)[0]
                           ] * len(table)
        final_table = final_table.append(table)
    final_table = final_table.reset_index()
    columns = list(final_table.columns)
    columns.remove("sample")
    columns = ["sample"] + columns
    final_table = final_table[columns]
    final_table.to_csv(outfile, index=False, sep="\t")


if __name__ == '__main__':
    run()
