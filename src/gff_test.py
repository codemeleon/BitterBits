from gengff import gff
# import pandas as pd

table, sequence =  gff.dataframe("/home/devil/tmp/112844.gff")
print(table.head())
