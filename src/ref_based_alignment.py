#!/usr/bin/env python
"""Feel free to use."""

import click
from os import system, path, makedirs
from Bio import SeqIO


@click.command()
@click.option("-ref", help="Reference sequence", type=str, default=None,
              show_default=True)
@click.option("-index", help="Generate index file", type=bool, default=False,
              show_default=True)
@click.option("-f", help="Forward sequences", type=str, default=None,
              show_default=True)
@click.option("-r", help="Reverse sequences", type=str, default=None,
              show_default=True)
@click.option("-ofile", help="Ourput file", type=str, default=None,
              show_default=True)
@click.option("-alnr", help="Alinger used",
              type=click.Choice(["bowtie2", "bwa"]),
              default="bowtie2", show_default=True)
@click.option("-clean", help="Delete itermediate file after finish task",
              type=bool, default=True, show_default=True)
def run(ref, index, f, r, ofile, alnr, clean):
    """Generate alternative sequence based given ref and reads."""
    assert ref, "Reference file is not given. Exiting ..."
    assert f, "Forward file not given. exiting ..."
    assert ofile, "Output file not given. Exiting ...."
    if not path.exists(ref):
        click.echo("Reference file path doesn't exist. Exiting")
        exit()
    if not path.exists(f):
        click.echo("Forward sequence file path doesn't exist. Exiting")
        exit()
    if not path.exists(r):
        click.echo("Reverse sequence file path doesn't exist. Exiting")
        exit()
    base_name = path.split(f)[1].split(".")[0]
    tmpf = "/tmp/%s" % base_name
    makedirs(tmpf)
    if alnr == "bowtie2":
        system("bowtie2-build -f %s %s/index" % (ref, tmpf))
        if r:
            system("bowtie2 -x %s/index -1 %s -2 %s -S %s/mapped.sam"
                   % (tmpf, f, r, tmpf))
        else:
            system("bowtie2 -x %s/index -U %s -S %s/mapped.sam"
                   % (tmpf, f, tmpf))

    system("samtools view -bS -o %s/mapped.bam %s/mapped.sam" %
           (tmpf, tmpf))
    system("samtools sort %s/mapped.bam %s/mapped_sorted" %
           (tmpf, tmpf))
    system("samtools index %s/mapped_sorted.bam" % (tmpf))
    system("""samtools mpileup -uf %s %s/mapped_sorted.bam |
           bcftools view -cg - |
           vcfutils.pl vcf2fq > %s/mapped.fq""" %
           (ref, tmpf, tmpf))
    ouput_file = open(ofile, "w")
    for rec in SeqIO.parse("%s/mapped.fq" % tmpf, "fastq"):
        ouput_file.write(">%s\n%s\n" % (rec.id, rec.seq))
    ouput_file.close()
    system("rm -rf %s" % tmpf)


if __name__ == '__main__':
    run()
