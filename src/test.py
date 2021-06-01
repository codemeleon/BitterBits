import matplotlib
import numpy as np
from pylab import *
from scipy.stats import pearsonr


2+2  # +"anmol"
matplotlib.use("Agg")

infile = open("/home/anjali/test/counts/rapamycin/counts.csv")
te1l = []
te2l = []
infile.readline()
for line in infile:
      linesplit = line.split(',')
      te1 = np.log2(max([float(linesplit[5]), 1])))
          te2=np.log2(max([float(linesplit[6]), 1])))
          te1l.append(te1)
          te2l.append(te2)
          infile.close()


          fig=figure()
          plt.scatter(te1l, te2l, s=0.2, color='grey')
          ylabel("Rapamycin(TE)#Replicate1", fontsize=10,
                 color='black', style='italic', fontweight='bold')
          xlabel("Rapamycin(TE)#Replicate2", fontsize=10,
                 color='black', style='italic', fontweight='bold')

          # scatter(gene_expression, te_var, s = 0.5,color = "#2E9AFE")


          coeff=round(pearsonr(te1l, te2l)[0], 3)
          text(2, 16, "Pearson coefficient=%s" %
               coeff, color='black', style='oblique')
          xticks([2, 4, 6, 8, 10, 12, 14, 16, 18],
                 [2, 4, 6, 8, 10, 12, 14, 16, 18])
          savefig("TE_rapa.png")
          clf()



          {'Abobo': "Ivory Coast",
           'AddisAbaba': 'Ethiopia',
           'Adjame': "Morocco",
           'Agadir': "Morocco",
           'Agnibilekrou',
           'Algeria',
           'Algiers',
           'Analavory',
           'Anjeva',
           'Annaba',
           'Antananarivo',
           'Antsirabe',
           'Attecoube',
           'Anmoehoririka',
           'AnmoeniMellal',
           'Anmoenimellal',
           'Anmolida',
           'Anmoouafle',
           'Anmooumerdes',
           'AnmourkinaFaso',
           'Anmoairo',
           'Anmoameroon',
           'CapeTown',
           'Capetown',
           'Casablanca',
           'CentralAfricanRepublic',
           'Cocody',
           'Congo',
           "CoteD'Ivoire",
           "Coted'Ivoire",
           "Coted'Ivorie",
           'DAKAR',
           'Dakar',
           'Djibouti',
           'Douala',
           'Durban',
           'Egypt',
           'Fes',
           'Fianarantsoa',
           'Gambia',
           'Ghana',
           'Hevecam',
           'Itaosy',
           'JOHANNESBURG',
           'Johanneburg',
           'Johannesburg',
           'KENYA',
           'Kenya',
           'Lusaka',
           'MOROCCO',
           'Madagascar',
           'Maevatanana',
           'Mahajanga',
           'Mali',
           'Mali071Ci',
           'Mali107Ci',
           'Manjakaray',
           'Marrakech',
           'Mauritius',
           'Meknes',
           'Moramanga',
           'Morocco',
           'Mozambique',
           'Niakhar',
           'Niger',
           'Nigeria',
           'NosyBe',
           'Nosybe',
           'Oujda',
           'Plateau',
           'Rabat',
           'Rwanda',
           'SOUTHAFRICA',
           'Sale',
           'Senegal',
           'Setif',
           'SierraLeone',
           'South-Africa',
           'SouthAfrica',
           'Southafrica',
           'Tanger',
           'Tanzania',
           'Temara',
           'Tiznit',
           'Toamasina',
           'Togo',
           'Tsaralalana',
           'Tsiroanomandidy',
           'Tunis',
           'Tunisia',
           'Uganda',
           'Yaounde',
           'Yopougon',
           'Zambia',
           'morocco',
           'southAfrica'}
atgtcgctgtttggagacacaattgcctacctgctttcattaacagaagatggagaaggcaaagcagaactagcagaaaaattacactgttggttcggtgggaaagaatttgacctagactctgccttggaatggataaaaaacaaaagatgcttaactgatatacagaaagcactaattggtgcctctatctgctttttaaaacccaaagaccaggaaagaaaaagaagattcatcacagagcccttatcgggaatgggaacaacagcaacaaaaaagaagggcctgattctggctgagagaaaaatgagaaaatgtgtgagctttcatgaagcatttgaaatagcagaaggccatgaaagctcagcgctactatattgtcttatggtcatgtacctgaatcctggaaattattcaatgcaagtaaaactagggacactctgtgctttgtgcgaaaaacaagcatcacattcacacagagctcatagcagagcagcgagatcctcagtgccaggagtgagacgggaaatgcagatggtctcagctatgaacacagcaaaaacaatgaatggaatgggaaaaggagaagacgtccaaaaactggcagaagaactgcaaagcaacattggagtattgagatcccttggggcaagtcagaagaacggggaaggaattgcaaaggatgtaatggaagtgctaaagcagagctctatgggaaattcagctcttgtgaagaaatacctataatgtttgaaccattccagattctttcgatttgttcttttattctatcagctctccatttcatggcttggacaataggacatttaaatcaaataaaaagaggagtaaacatgaaaataagaataaaggggccaaataaagagacaataaacagagaggtatcaattttgagacacagttaccaaaaagaaattcaggctaaggaagcaatgaaggaagtactctctgacaacatggaagtattgagtgaccacatagtaattgaggggctttctgctgaagagataataaaaatgggtgaaacagttttggaggtagaagaatcgcattaa

CTTGAGACTATCGAACTCTGTAAACAAAGACAGACAAGCTACCAATTNGCTCCTCTACTTATTCTTGTAGGAAATGCTGGAGCCAAATTTCTAGACAGAG
