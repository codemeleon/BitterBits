import pysam
import pandas as pd
from glob import glob
from multiprocessing import Pool
import numpy as np
from scipy import stats
from os import path



def run(fl):
    comments = ""
    with open(fl) as c:
        for line in c:
            if not line.startswith('#'):
                break
            comments += line
    comments = comments[:-1]+"\tendbases\n"
    
    data = pd.read_table(fl, comment='#', header=None)
    samfile = pysam.AlignmentFile(f"{fl.split('_fill')[0]}_L001001_sequence-alignment.bam", "rb")
    ends_bases_count = []
    for _, row in data.iterrows():
        for pileupcolumn in samfile.pileup("NC_045512.3", row[1]-2,row[1], min_base_quality=0,max_depth=50000):
            if pileupcolumn.pos == row[1]-1:
                n = []
                nc = 0
                for pileupread in pileupcolumn.pileups:
                    if not pileupread.is_del and not pileupread.is_refskip:
                        if (pileupread.alignment.query_sequence[pileupread.query_position] != row[3]):
                            nc += 1
                            if ((pileupread.query_position < 10) or ((pileupread.alignment.qlen - (pileupread.query_position+1))<10)):
#                             n += 1 # Later check for consistancy
                                n.append(min((pileupread.query_position, pileupread.alignment.qlen - (pileupread.query_position+1))))
#                 print(stats.mode(n).count[0])
                if not n:
                    most_common_count = 0
                else:
                    most_common_count = stats.mode(n).count[0]
#                 print(n, stats.mode(n).count[0],most_common_count, nc, most_common_count/nc)
                ends_bases_count.append(most_common_count/nc)
                break
    data["endsbases"] = ends_bases_count
    f = open(f"res/{path.split(fl)[1]}", 'w')
    f.write(comments)
    data.to_csv(f, header=False, index=False, sep="\t")
    f.close()
    
    # data.to_csv(f"res/{path.split(fl)[1]}", header=False, index=False)
                            
p = Pool(24)

_ = p.map(run, glob("ans/*_fill.vcf"))