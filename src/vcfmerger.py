import gzip
from glob import glob
from os import path

import click
import numpy as np
import pandas as pd


def depth(x):
    for y in x.split(";"):
        z = y.split("=")
        if z[0] == "DP4":
            return sum(list(map(int, z[1].split(','))))


@click.command()
@click.option(
    "-vcfd",
    help="VCF files containing folder",
    type=str,
    default="vcffiles",
    show_default=True)
@click.option(
    "-sfx",
    help="Gff file suffix",
    type=str,
    default=".vcf.gz",
    show_default=True)
@click.option(
    "-dpth",
    help="Minimum good reads depth for SNP",
    type=int,
    default=30,
    show_default=True)
@click.option(
    "-typ",
    help="Type of alteration to report",
    type=click.Choice(["snp", "indels", "both"]),
    default="both",
    show_default=True)
@click.option(
    "-opx",
    help="Output file prefix",
    type=str,
    default="output",
    show_default=True)
def run(vcfd, sfx, dpth, typ, opx):
    """
    This scrimt to merge vcf files base on the reported variant coordinates and
    the given read depth. It provies, genomic coordinates and alternative bases

    The File must contain following columns, considering chromosome is always
    the same.
    CHROM at col 1
    POS at col 2
    ALT at col 5
    INFO at col8
    INFO must contain DP4 variable (Required for depth calculation)
    """
    table_indel = pd.DataFrame.from_dict({
        'POS': []
    })  # Keep merging without waiting to avoid memory usages
    table_snp = pd.DataFrame.from_dict({
        'POS': []
    })  # Keep merging without waiting to avoid memory usages
    # TODO: Indel SNP coordinate overlap
    # fn_open = gzip.open if filename.endswith('.gz') else open
    # use "compression" in read_table
    # print("Running First Step")

    for i, fl in enumerate(glob("%s/*%s" % (vcfd, sfx))):
        fl_bs = path.split(fl)[1].split('.', 1)[0]
        vcffile = pd.read_table(
            fl,
            header=None,
            comment='#',
            compression="gzip",
            usecols=[0, 1, 4, 7])
        vcffile = vcffile.rename(columns={
            0: 'CHROM',
            1: 'POS',
            4: fl_bs,
            7: 'INFO'
        })
        vcffile["POS"] = list(map(str, vcffile["POS"]))
        vcffile["POS"] = vcffile[['CHROM', 'POS']].apply(
            lambda x: '.'.join(x), axis=1)
        vcffile["DEPTH"] = list(map(depth, vcffile["INFO"]))
        vcffile = vcffile[vcffile["DEPTH"] > dpth]
        vcffile["TYPE"] = "SNP"
        vcffile.loc[vcffile["INFO"].str.contains("INDEL"), "TYPE"] = "INDEL"

        vcf_snp = vcffile.loc[vcffile["TYPE"] == "SNP", ['POS', fl_bs]]
        vcf_indel = vcffile.loc[vcffile["TYPE"] == "INDEL", ['POS', fl_bs]]
        # print(vcf_snp)
        # break
        del vcffile  # ["INFO"]

        table_snp = pd.merge(table_snp, vcf_snp, on="POS", how='outer')
        table_indel = pd.merge(table_indel, vcf_indel, on="POS", how='outer')
        # print(i, fl, table_snp.shape, table_indel.shape)
    # print("Running second step")
    if typ in ['snp', 'both']:
        table_snp.to_csv("%s.snp.tsv" % opx, index=False, sep="\t")
    # print("Running third step")
    if typ in ['indel', 'both']:
        table_indel.to_csv("%s.indel.tsv" % opx, index=False, sep="\t")

    # print("Running fourth step")
    if typ == 'both':
        common = list(set(table_snp["POS"]) & set(table_indel["POS"]))
        common.sort()
        table_indel_c = table_indel[table_indel["POS"].isin(common)]
        table_snp_c = table_snp[table_snp["POS"].isin(common)]
        with open("%s.snp.indel.common.coor" % opx, "w") as fout:
            fout.write("POS\tSNP\tINDEL\n")
            for c in common:
                snp = table_snp_c.loc[table_snp_c["POS"] == c, ]
                snp = snp.columns[~snp.isnull().any()]
                snp = set(snp) - set(["POS"])
                snp = ",".join(snp)
                indel = table_indel_c.loc[table_indel_c["POS"] == c, ]
                indel = indel.columns[~indel.isnull().any()]
                indel = set(indel) - set(["POS"])
                indel = ",".join(indel)

                fout.write("%s\t%s\t%s\n" % (c, snp, indel))
        # for _, row in table_indel.iterrows():
        #     print(row[~pd.isnull(row)].index)
        # print(row.index[pd.isnull()])

    # print(table_indel.head())


if __name__ == '__main__':
    run()
