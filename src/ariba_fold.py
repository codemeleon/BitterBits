#!/usr/bin/env python
from sys import argv
from os import path, makedirs
from glob import glob
from os import system, path

if path.exists(argv[3]):
    if path.isfile(argv[3]):
        print("Given Output Path is a file. Exiting ...")
    if not path.isdir(argv[3]):
        makedirs(argv[3])

system("ariba prepareref --all_coding yes -f %s tmp_ariba" % argv[2])

for fl in glob("%s/*_1.fastq.gz" % argv[1]):
    file_base = path.split(fl)[1].split("_1.fa")[0]
    path_base = fl.split("_1.fa")[0]

    # print(file_base)
    system("bsub -e %s/%s/%s.e -o %s/%s/%s.o -R 'select[mem>4000] rusage[mem=4000]"
           " span[hosts=1]' -n 1 -M4000 ariba run tmp_ariba %s_1.fastq.gz"
           " %s_2.fastq.gz %s/out.%s"
           % (argv[3], file_base, file_base, argv[3], file_base, file_base,
              path_base, path_base, argv[3], file_base))
