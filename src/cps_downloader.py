from Bio import Entrez
import pandas as pd
import time
# Might not work completely
Entrez.email='manu_3601@yahoo.com'



sero_table = pd.read_csv("data/CPS_strain_list.txt")

for _, row in sero_table.iterrows():
    # print(row)
    handle = Entrez.esearch(db="nucleotide", retmax=10,
                            term='''"%s"[All Fields] AND
                            "Streptococcus pneumoniae"[porgn] AND
                            bacteria[filter]''' % row["Strain"], idtype="acc")
    record = Entrez.read(handle)
    time.sleep(3)
    if len(record["IdList"]) != 1:
        print(row)
        continue
    handle = Entrez.efetch(db="nucleotide", id=record["IdList"][0], rettype="gb", retmode="text")
    with open("%s.gb" % row["Serotype"], "w") as fout:
        fout.write(handle.read())
    time.sleep(3)
