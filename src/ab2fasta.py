from Bio import SeqIO
import click
from os import path
from glob import glob


def ab1_to_fasta(ab1):
    """Extract numceotide sequences from ab1 file
    :ab1: ab1 sanger trace file
    :returns: Nucleotide sequence
    """
    record = SeqIO.read(ab1, "abi")
    return "".join([chr(c) for c in record.annotations["abif_raw"]["PBAS1"]])


@click.command()
@click.option(
    "-ab1",
    help="ab1 file or folder containing ab1 files",
    type=str,
    default=None,
    show_default=True)
@click.option(
    "-res",
    help="Output folder",
    type=str,
    default=".",
    show_default=True)
def run(ab1, res):
    """Generate ab1 equivalent fasta

    :ab1: fasta or folder
    :returns: None

    """
    if not path:
        exit("Input path not given")
    if not path.exists(ab1):
        exit(f"Given {ab1} path doesn't exist")
    if path.isdir(ab1):
        for fl in glob(f"{ab1}/*.ab1"):
            seq = ab1_to_fasta(fl)
            flb = path.split(fl)[1].split(".ab1")[0]
            flbf = f"{res}/{flb}.fasta"
            with open(flbf, "w") as fout:
                fout.write(f">{flb}\n{seq}\n")
    elif path.isfile(ab1):
        seq = ab1_to_fasta(ab1)
        flb = path.split(ab1)[1].split(".ab1")[0]
        flbf = f"{res}/{flb}.fasta"
        with open(flb, "w") as fout:
            fout.write(f">{flb}\n{seq}\n")


if __name__ == "__main__":
    run()
