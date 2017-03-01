#!/usr/bin/env python
import click
import os
from subprocess import call, Popen, PIPE, STDOUT
from os import remove, path
import pandas as pd
from Bio import SeqIO

# http://click.pocoo.org/5/advanced/

# Write program to check the existance of program or program path


def is_tool(name):
    """Repopring if the tool exist in the path."""
    # Code taken from
    # "http://stackoverflow.com/questions/11210104/check-if-a-program-exists-from-a-python-script"
    try:
        devnull = open(os.devnull)
        Popen([name, '-h'], stdout=devnull, stderr=devnull).communicate()
    except OSError as err:
        if err.errno == os.errno.ENOENT:
            return False
    return True


def find_duplicates(fl):
    """Considering sequences are in all upper or lower case."""
    'File is in fasta format'
    seq_dict = SeqIO.to_dict(SeqIO.parse(fl, 'fasta'))
    df = pd.DataFrame.from_dict({'name': list(seq_dict.keys()),
                                 'seq': list(seq_dict.values())})
    df['seq'] = df['seq'].apply(''.join)
    # Pirated from
    "http://stackoverflow.com/questions/19530568/can-pandas-groupby-aggregate\
    guilsinglright-into-a-list-rather-than-sum-mean-etc"
    df = df.groupby('seq')
    df = df.aggregate(lambda x: tuple(x))
    # -----------------------------------------------------------------------
    return df


def aligner(infile=None, outfile=None, alntool='clustalw', toolpath=None,
            additional_arg=""):
    if not infile:
        click.echo("Input file not given!!")
        exit(1)
    elif not outfile:
        click.echo("Output file not given")
        exit(1)
    elif not alntool:
        click.echo("Alignment tool not given")
        exit(1)
    elif alntool not in ['clustalw', 'muscle', 'prank']:
        click.echo("Alignment tool is not from the list ['clustalw', 'muscle',\
                   'prank']")
        exit(1)
    else:
        if not toolpath:
            toolpath = alntool
        if not path.isfile(infile):
            click.echo("Input file doesn't exist")
            exit(1)
        if not is_tool(toolpath):
            click.echo("Given tool path or tool in path doesn't exist")
            exit(1)

        # Check the program

    df = find_duplicates(infile)
    if df.shape[0] == 1:
        with open(outfile, "w") as fout:
            fout.write(open(infile).read())
        return
    tmp_in = "%s.tmp_in" % infile
    tmp_out = "%s.tmp_out" % infile
    seq_groups = {}
    with open(tmp_in, "w") as ofile:
        for i, row in df.iterrows():
            seq_groups[row['name'][0]] = row['name']
            ofile.write(">%s\n%s\n" % (row['name'][0], i))
    # Check if the tool exists
    # toolpath = ""  # alter it. For now use it for avoiding the error
    # Only fasta format please (in and out)
    if additional_arg != "":
        command = {
            'muscle': [toolpath, additional_arg, '-in', tmp_in,
                       '-out', tmp_out],
            'clustalw': [toolpath, additional_arg, '-align', '-output=fasta',
                         '-infile=%s' % tmp_in, '-outfile=%s' % tmp_out],
            'prank': [toolpath, additional_arg, '-d=%s' % tmp_in,
                      '-o=%s' % tmp_out]
        }
    else:
        command = {
            'muscle': [toolpath, '-in', tmp_in, '-out', tmp_out],
            'clustalw': [toolpath, '-align', '-output=fasta',
                         '-infile=%s' % tmp_in, '-outfile=%s' % tmp_out],
            'prank': [toolpath, '-d=%s' % tmp_in,
                      '-o=%s' % tmp_out]
        }
    Popen(command[alntool], stdout=PIPE, stderr=STDOUT).communicate()

    with open(outfile, "w") as fout:
        for rec in SeqIO.parse(tmp_out, 'fasta'):
            for gid in seq_groups[rec.id]:
                fout.write(">%s\n%s\n" % (gid, rec.seq))
    remove(tmp_in)
    remove(tmp_out)


@click.command()
@click.option("--i", help="Input fasta file", type=str,
              default=None, show_default=True)
@click.option("--o", help="Output fasta file. please don't change to other\
              format", default=None, show_default=True)
@click.option("--alntool", help="Alignment tool",
              type=click.Choice(['clustalw', 'muscle', 'prank']),
              default='muscle', show_default=True)
@click.option("--toolpath", help="Selected toolpath. If path not given, tool\
              will be considered in system path", default=None,
              show_default=True)
@click.option("--additional_arg", help="Additional arguments for the " +
              "selected alignment tool. Please check help page for the " +
              "selected program Please keepin mind, output must be fasta",
              type=str, default="", show_default=True)
def run(i, o, alntool, toolpath, additional_arg):
    r"""Program boosts the speed of alignments where sequences exists \
    in multiple copies."""
    aligner(i, o, alntool, toolpath, additional_arg)


if __name__ == '__main__':
    run()
