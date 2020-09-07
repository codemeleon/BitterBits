#!/usr/bin/env python

# This is a python wrapper for virulantfinder
# Install the virulantfinder docker from
# https://bitbucket.org/genomicepidemiology/virulencefinder/src/master/

# Download and contruct databases by flowing the Instrusctions at above page


# Before execution of this program, make docker available for the user
# Please follow steps 1-3 at https://www.thegeekdiary.com/run-docker-as-a-non-root-user/

from glob import glob
from os import makedirs, path, system

import click


@click.command()
@click.option(
    "-od",
    help="Output directory",
    type=str,
    default="vf_result",
    show_default=True,
)
@click.option(
    "-id", help="Input directory", type=str, default="..", show_default=False
)
@click.option(
    "-ft",
    help="Filetype",
    type=click.Choice(["fasta", "gff"], case_sensitive=True),
    default="gff",
    show_default=True,
)  # Add list
#  @click.option()
def run(od, id, ft):
    if not id:
        exit("Input folder is not given. Exiting . . . . .")
    if not path.isdir(id):
        exit("Given input folder path is not a folder. Exiting . . . .")
    if not od:
        exit("Outputfile is not given. Exit . . . .")
    # Add secondory command as specie
    for fl in glob(f"{id}/*.{ft}"):
        print(fl)
        flb = path.split(fl)[1].split(".")[0]
        if ft == "gff":
            seq = open(fl).read().split("##FASTA")[1]
            with open("tmp.fasta", "w") as fout:
                fout.write(seq.strip())

        else:
            with open("tmp.fasta", "w") as fout:
                fout.write(open(fl).read())
        if not path.exists(f"{od}/{flb}"):
            makedirs(f"{od}/{flb}", exist_ok=True)

        outd = f"{od}/{flb}"

        system(
            f"docker run --rm -it -v "
            "/.anmol/databases/virulencefinder_db:/database "
            f"-v {path.abspath('.')}:/workdir "
            f"virulencefinder -i {path.relpath('tmp.fasta')} "
            f"-o {path.relpath(outd)} -q"
        )

        #  break


if __name__ == "__main__":
    run()
