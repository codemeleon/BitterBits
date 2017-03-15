import sys
from Bio import Entrez
import time
Entrez.email = "manu_3601@yahoo.com"
# Modify the code to make it more usable
id_1 = sys.stdin.readline().strip()
while id_1:
    record = Entrez.read(Entrez.elink(dbfrom="protein", db="nuccore", id=id_1))

    id_2 = record[0]['LinkSetDb'][0]['Link'][0]["Id"]
    time.sleep(4)
    with open("%s_%s.gb" % (id_1, id_2), "w") as fout:
        record = Entrez.efetch(db="nucleotide", id=id_2,
                               rettype="gbwithparts", retmode="text").read()
        fout.write(record)
    time.sleep(4)
    id_1 = sys.stdin.readline().strip()
