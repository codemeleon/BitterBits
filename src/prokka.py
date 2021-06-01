#!/usr/bin/env python
"""This script to annotated all fasta files in the folder considering all are
bacterial. Puts All the results in Prokka folder"""


from glob import glob
from os import makedirs, path
from subprocess import Popen

from multiprocess import cpu_count

cpu_use = int(cpu_count() * 0.75)

#  makedirs("Prokka")


for fl in glob("*"):
    flb = fl.split(".")[0]
    pkk_cmd = [
        "prokka",
        fl,
        "--cpus",
        f"{cpu_use}",
        "--outdir",
        f"Prokka/{flb}",
        "--prefix",
        flb,
        "--locustag",
        flb,
    ]
    Popen(pkk_cmd).communicate()
