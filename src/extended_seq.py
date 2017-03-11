#!/usr/bin/env python
"""Generate sequences for given ids. Requires gengff folder."""
import click
from gengff import gff


@click.command()
@click.option("-seqids", help="Comma seperated ids", type=str, default=None,
              show_default=True)
@click.option("-nttfl", help="Annotation file including dna sequences",
              type=str, default=None, show_default=True)
@click.option("-nttfmt", help="Annotation file format",
              type=click.Choice(["gff", "genbank"]), default="gff",
              show_default=True)
@click.option("-flnk", help="Flanking region", type=int, default=1000,
              show_default=True)
@click.option("-mindist", help="Minimum nt distance to consider them together",
              type=int, defaut=500, show_default=True)
@click.option("-outfl", help="Outfile. default:stdout", type=str, default=None,
              show_default=True)
def run(seqids, nttfl, nttfmt, flnk, outfl, mindist):
    """Generate sequences for given ids."""
    assert seqids, "Seqids are not given. Exiting ...."
    assert nttfl, "Annotaton file not given. Exiting ...."
    seqids = seqids.split(",")
    if nttfmt in ["gff", "gtf"]:
        table, sequences = gff.dataframe(nttfl)
    table = table[table["locus_tag"].isin(seqids)
                  ][["locus_tag", "start", "end", "seqname"]]
    outfile = open(outfl, "w")
    for seqname in set(table["seqname"]):
        sub_table = table[table["seqname"] == seqname]
        if len(sub_table) == 1:
            for i, row in sub_table.iterrows():
                outfile.write(">%s\n%s\n" %
                              (row["locus_tag"],
                               sequences[seqname][row["start"] - flnk:
                                                  row["end"] + flnk]))
        else:
            sub_table = sub_table.sort_values("start").reset_index()
            start = 0
            tag_name = []
            for i, row in sub_table.iterrows():
                if i != len(sub_table) - 1:
                    if abs(sub_table.ix[i+1]["start"] - row["end"]) < mindist:
                        if not start:
                            start = row["start"]
                        tag_name.append(row["locus_tag"])
                    else:
                        if not start:
                            outfile.write(">%s\n%s\n" %
                                          (row["locus_tag"],
                                           sequences[seqname
                                                     ][row["start"] - flnk:
                                                       row["end"] + flnk]))
                        else:
                            outfile.write(">%s\n%s\n" %
                                          ('___'.join(tag_name),
                                           sequences[seqname
                                                     ][start - flnk:
                                                       row["end"] + flnk]))
                        start = 0
                        tag_name = []
                else:
                    if not start:
                        outfile.write(">%s\n%s\n" %
                                      (row["locus_tag"],
                                       sequences[seqname
                                                 ][row["start"] - flnk:
                                                   row["end"] + flnk]))
                    else:
                        outfile.write(">%s\n%s\n" %
                                      ('___'.join(tag_name),
                                       sequences[seqname
                                                 ][start - flnk:
                                                   row["end"] + flnk]))
    outfile.close()


if __name__ == '__main__':
    run()
