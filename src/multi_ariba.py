#!/nfs/users/nfs_a/ak20/miniconda3/bin/python
from glob import glob
import click
from os import system, path, makedirs


@click.command()
@click.option(
    "-dbf",
    help="Preprocessed reference database",
    type=str,
    default=None,
    show_default=True)
@click.option(
    "-fqf", help="Fastq folder", type=str, default=None, show_default=True)
@click.option(
    "-fr",
    help="Forward mate representation",
    type=str,
    default="_1",
    show_default=True)
@click.option(
    "-rv",
    help="Reverse mate representations",
    type=str,
    default="_2",
    show_default=True)
@click.option(
    "-sfx",
    help="Fastq file sufix",
    type=str,
    default=".fastq.gz",
    show_default=True)
@click.option(
    "-bsb",
    help="Bsub separate for each file",
    type=bool,
    default=True,
    show_default=True)
@click.option(
    "-outd", help="Outpur directory", type=str, default=".", show_default=True)
def run(dbf, fqf, fr, rv, sfx, bsb, outd):
    """This program is to run ariba over whole directory files."""
    if not fqf:
        exit("Fastq folder not given. Exiting . . . .")
    if not path.isdir(fqf):
        exit("Given Fastq folder path is not directory. Exiting . . .")
    if outd != ".":
        makedirs(outd, exist_ok=True)
    if bsb:
        for f1 in glob("%s/*%s%s" % (fqf, fr, sfx)):
            bs = path.split(f1)[1].split("%s%s" % (fr, sfx))[0]
            f2 = f1.replace("%s%s" % (fr, sfx), "%s%s" % (rv, sfx))
            system("~sh16/scripts/bsubber.py -o log.o \
                    -e log.e -m 5 \
                    /nfs/users/nfs_a/ak20/miniconda3/bin/ariba run %s \
                    %s %s %s/%s" % (dbf, f1, f2, outd, bs))
    else:
        for f1 in glob("%s/*%s%s" % (fqf, fr, sfx)):
            bs = path.split(f1)[1].split("%s%s" % (fr, sfx))[0]
            f2 = f1.replace("%s%s" % (fr, sfx), "%s%s" % (rv, sfx))
            system("/nfs/users/nfs_a/ak20/miniconda3/bin/ariba run %s \
                    %s %s %s/%s" % (dbf, f1, f2, outd, bs))


if __name__ == '__main__':
    run()
