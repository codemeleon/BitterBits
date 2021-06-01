import os
import urllib
from time import sleep

import click
from Bio import Entrez

# This code is taken from https://dmnfarrell.github.io/bioinformatics/assemblies-genbank-python
# and modified for my own need

def get_assembly_summary(id):
    """Get esummary for an entrez id"""
    esummary_handle = Entrez.esummary(db="assembly", id=id, report="full")
    esummary_record = Entrez.read(esummary_handle)
    return esummary_record


@click.command()
@click.option("-e", help="Email id", type=str, default=None, show_default=True)
@click.option("-f", help="File having information one per line", type=str, default=None, show_default=True)
def run(e, f):
    if not e:
        exit("Email is not given")
    if not f:
        exit("Input file not given")

    with open(f) as fin:
        for line in fin:
            if line == '\n':
                continue
            summary = get_assembly_summary(line[:-1])

            url = summary['DocumentSummarySet']['DocumentSummary'][0]['FtpPath_RefSeq']
            if url == '':
                continue
            label = os.path.basename(url)
            link = os.path.join(url,label+'_genomic.fna.gz')
            urllib.request.urlretrieve(link, f'{line[:-1]}.fna.gz')
            sleep(10)



if __name__ == '__main__':
    run()
