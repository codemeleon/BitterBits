from glob import glob
from sys import argv


for fl in glob(argv[2]):
    with open(fl) as fin:
        bl = 0
        for line in fin:
            if argv[1] in line:
                bl = 1
                break
        if not bl:
            print(fl)