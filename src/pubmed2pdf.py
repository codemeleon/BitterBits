#!/usr/bin/env python

import urllib3
import requests
from bs4 import BeautifulSoup as bs
import click
from os import path, makedirs
from glob import glob
from time import sleep

headers = requests.utils.default_headers()
headers['User-Agent'] = 'Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/56.0.2924.87 Safari/537.36'


def pdfload(outd, pid):
    http = urllib3.PoolManager()
    page = http.request('GET', f"https://pubmed.ncbi.nlm.nih.gov/{pid}/")
    print(f"Featching {pid} data")
    soup = bs(page.data, "html.parser")
    lks = soup.find_all("a", attrs={"data-ga-category":"full_text"})
    links = {}
    print(f"Extracting PMCID and DOI for {pid}")
    print(lks)
    # return
    for link in lks:
         if link.attrs['data-ga-action'] == "PMID":
             continue
         links[link.attrs['data-ga-action']] = link.attrs['href']
    if "PMCID" in links:
        pmc_data = http.request('GET', links['PMCID'])
        pmc_soup = bs(pmc_data.data, "html.parser")
        plinks = []
        print(plinks)
        try:
            for link in pmc_soup.find_all("a") + pmc_soup.find_all("link"):
                if link.attrs['href'].lower().endswith(".pdf"):
                    plinks.append("https://www.ncbi.nlm.nih.gov"+link.attrs['href'])
            plinks = list(set(plinks))
            if len(plinks) ==0 :
                return pid
            if len(plinks) == 1:
                r = requests.get(plinks[0], headers=headers)
                with open(f"{outd}/{pid}.pdf", 'wb') as f:
                    f.write(r.content)
            else:
                for i, lk in enumerate(plinks):
                    r = requests.get(lk, headers=headers)
                    with open(f"{outd}/{pid}_{i}.pdf", 'wb') as f:
                        f.write(r.content)
                    sleep(5)
        except:
            return pid
    elif "DOI" in links:
        try:
            alt_data = http.request('GET', links['DOI'],headers=headers)
            alt_soup = BeautifulSoup(alt_data.data, "html.parser")
            alt_links = []
            for link in alt_soup.find_all("a"):
                if link.attrs['href'].lower().endswith(".pdf"):
                    alt_links.append(link.attrs['href'])
            plinks = list(set(alt_links))
            if len(plinks) ==0 :
                return pid
            if len(plinks) == 1:
                r = requests.get(plinks[0], headers=headers)
                with open(f"{outd}/{pid}.pdf", 'wb') as f:
                    f.write(r.content)
            else:
                for i, lk in enumerate(plinks):
                    r = requests.get(lk, headers=headers)
                    with open(f"{outd}/{pid}_{i}.pdf", 'wb') as f:
                        f.write(r.content)
                    sleep(5)
        except:
            return pid




@click.command()
@click.option("-pid", help="Pubmed ID, comma seprated", type=str, default=None, show_default=True)
@click.option("-pmf", help="Pubmed ID file, one id per line", type=str, default=None, show_default=True)
@click.option("-outd", help="Output Directory", type=str, default="fetched_pdfs", show_default=True)
@click.option("-uid", help="Unfetch id file", type=str, default="unfetched_ids", show_default=True)
def main(pid, pmf, outd, uid):
    """
    This script allows to download pdf copy of papers. First it search PMC.
    If the pdf copy is not available there. Then it tries to  journal websites.
    This work is still in progress.
    """
    if not path.isdir(outd):
        makedirs(outd, exist_ok=True)
    if not pid and not pmf:
        exit("ID or ID file is not given. Exiting .......")
    elif pid and pmf:
        exit("ID and ID file,  both not given. Use only one. Exiting .......")
    elif pid:
        pids = [p.strip() for p in pid.split(",")]
    else:
        if not path.isfile(pmf):
            exit("Given file doesn't exist. Exiting ....")
        else:
            pids = [p.strip() for p in open(pmf).readlines()]
            pids = [p for p in pids if len(p)>4]
    uf = open(uid,"w")
    for p in pids:
        if len(glob(f"{outd}/{p}_*.pdf")) or path.isfile(f"{outd}/{p}.pdf"):
            print(f"{p} already fetched. Skipping .....")
            continue
        notfound = pdfload(outd, p)
        if notfound:
            uf.write(f"{notfound}\n")
        sleep(5)
    uf.close()




if __name__ == '__main__':
    main()
