import click
from os import path, makedirs
from bs4 import BeautifulSoup
import urllib.request as ur
import requests

page = requests.get("http://www.ebi.ac.uk/ena/data/view/ERR018092&display=html")
soup = BeautifulSoup(page)
help(page)
page.text()
soup.unicode

table1 = soup.find("table")
table1

@click.command()
@click.option("-f", help="ID file. One each line", type=str,
              default=None, show_default=True)
@click.option("-d", help="Output directory", type=str,
              default=None, show_default=True)
def run(f, d):
    """Download Fastq sequences from ENA."""
    if not path.exists(f) or not path.isfile(f):
        click.echo("Given file or filepath doesn't exist.")
        exit(1)
    if not path.isdir(d):
        click.echo("Given folder doesn't exist. Creating ...")
        makedirs(d)

    pass


if __name__ == '__main__':
    run()
NTTCCAACCAATCAATAAAATTACGTTCAATCGTCACAGCCAGATCATAAAGTGCAGAGTTTCGGTCAGATAAACCAAAATCAATTACGGCAGTAATTTCAGCTTTGGCATCCGGCGCTGACCAAAACAGATTAGAAGCGTGTAAATCGT
