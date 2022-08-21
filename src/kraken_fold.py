from glob import glob
from multiprocessing import Pool
from os import path
from subprocess import Popen


def kraken(samp):
    """TODO: Docstring for kraken.

    :fin: TODO
    :returns: TODO

    """
    kdb = "/.anmol/databases/minikraken_20171019_8GB"
    fwd = "/data/Data/ESKAPE_Malawi/Resistance_Genes/fastq_files/{}_1.fastq.gz".format(
        samp
    )
    rev = "/data/Data/ESKAPE_Malawi/Resistance_Genes/fastq_files/{}_2.fastq.gz".format(
        samp
    )

    if not path.exists(fwd):
        print("{} does not exist".format(fwd))
        return
    if not path.exists(rev):
        print("{} does not exist".format(rev))
        return
    kraken_cmd1 = [
        "kraken",
        "--threads",
        "1",
        "--fastq-input",
        "--gzip-compressed",
        "--db",
        kdb,
        "--output",
        samp,
        "--paired",
        fwd,
        rev,
    ]
    print(" ".join(kraken_cmd1))
    kraken_cmd2 = ["kraken-report", "--db", kdb, samp]
    Popen(kraken_cmd1).communicate()
    with open(f"{samp}.kraken", "w") as of:
        Popen(kraken_cmd2, stdout=of).communicate()


# files = glob("/data/Data/ESKAPE_Malawi/Resistance_Genes/*_fastq.gz")
ids = None

with open("samples") as samples:
    ids = [line.strip() for line in samples]


p = Pool(20)
p.map(kraken, ids[:2])
