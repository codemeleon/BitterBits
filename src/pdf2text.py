#!/usr/bin/env python
"""This script extract text from all the pdf file in current folder name
'Library' and generated markdown files."""

# Please install pdfminer.six using conda or pip

from glob import glob
from os import path, system


def main():
    print(path.abspath("."))
    if path.split(path.abspath("."))[1] != "Library":
        exit("Not a 'Library' folder")

    for fl in glob("*.pdf"):
        print(f"Parsing {fl}")
        outfile = fl.replace(".pdf", ".md")
        if path.exists(outfile):
            continue
        else:
            system(f"pdf2txt.py {fl} |"
                   "awk '{if (length($0)==0 || (length($0)>10)) print $0}' |"
                   "grep -A1 . |"
                   "grep -v '^--$' >"
                   f"{outfile}")


if __name__ == '__main__':
    main()
