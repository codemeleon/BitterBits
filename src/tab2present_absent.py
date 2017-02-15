#!/usr/bin/env python

"""Present Absent table"""

import pandas as pd
import click


@click.command()
@click.option("-intab", help="Input Table",
              type=str, default=None, show_default=True)
@click.option("-outtab", help="Output table",
              type=str, deafult=None, show_default=True)
@click.option("-min_samp", help="Minimun sample count",
              type=int, default=1, show_deafult=True)
@click.option("-max_samp", help="Max number of samples [-1 for all]",
              type=int, default=-1, show_default=True)
@click.option("-rot", help="Ratate table", type=bool,
              deafult=True, show_deafult=True)
def run(intab, outtab, min_samp, max_samp, rot):
    """Generate present absent table."""
    assert intab, "Input table file not give. Exiting ...."
    assert outtab, "Output table file not give. Exiting ...."
    assert pd.read_table(intab), "Table file doesn't exist or not fomatted"
    in_tab = pd.read_table(intab, index_col="cluster")
    if max_samp > 0:
        in_tab = in_tab[in_tab["samp_count"] < max_samp]
    in_tab = in_tab[in_tab["samp_count"] < min_samp]
    for k in ["samp_count", "seq_count", "min",
              "mean", "sd", "median", "max"]:
        del in_tab[k]
    in_tab[in_tab != "*"] = 1
    in_tab[in_tab == "*"] = 0
    if rot:
        in_tab = in_tab.T
    in_tab.to_csv(outtab, sep="\t")


if  __name__ == '__main__':
    run()
