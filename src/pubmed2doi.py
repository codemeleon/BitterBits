#!/usr/bin/env python

import urllib3
import requests
from bs4 import BeautifulSoup as bs
import click
from time import sleep

headers = requests.utils.default_headers()
headers[
    "User-Agent"
] = "Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/56.0.2924.87 Safari/537.36"


def pdfload(pid):
    http = urllib3.PoolManager()
    page = http.request("GET", f"https://pubmed.ncbi.nlm.nih.gov/{pid}/")
    soup = bs(page.data, "html.parser")
    lks = soup.find_all("a", attrs={"data-ga-category": "full_text"})
    links = {}
    for link in lks:
        links[link.attrs["data-ga-action"]] = link.attrs["href"]
    if "DOI" in links:
        print(links["DOI"].split("doi.org/")[1])
        sleep(5)
    else:
        print("DOI not found")


@click.command()
@click.option(
    "-pid",
    help="Pubmed ID, comma seprated",
    type=str,
    default="28079170",
    show_default=True,
)
def main(pid):
    """
    This script allows to download pdf copy of papers. First it search PMC.
    If the pdf copy is not available there. Then it tries to  journal websites.
    This work is still in progress.
    """
    if not pid:
        exit("ID is not given. Exiting .......")
    pids = [p.strip() for p in pid.split(",")]
    for p in pids:
        pdfload(p)


if __name__ == "__main__":
    main()
