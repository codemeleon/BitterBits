#!/usr/bin/env python
import click
from glob import glob
from os import path, system

@click.command()
@click.option("-f1", help="First gff Folder", type=str, default=None,
              show_default=True)
@click.option("-f2", help="Second gff Folder", type=str, default=None,
              show_default=True)
@click.option("-out", help="Output file", type=str, default=None,
              show_default=True)
def run(f1, f2, out):
    if not path.isdir(f1):
        click.echo("First folder doesn't exists. Exiting ..")
        exit(1)
    if not path.isdir(f2):
        click.echo("Secons folder doesn't exists. Exiting ..")
        exit(1)
    if not out:
        click.echo("Outfile not given. Exiting ..")
        exit(1)
    fl1 = ','.join(glob("%s/*.gff" %f1))
    fl2 = ','.join(glob("%s/*.gff" %f2))
    system("query_pan_genome -a difference --input_set_one %s --input_set_two"
           " %s -o %s" %(fl1, fl2, out))

if __name__ == '__main__':
    run()
