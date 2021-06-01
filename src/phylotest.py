from Bio import Phylo

tree = Phylo.read("/home/devil/Documents/Influenza/"
                  "NewAnalysis/Trees/RAxML_bestTree"
                  ".FluB_iva_HA_aln_trim_select", "newick")
tree.ladderize()
tree.rooted()
Phylo.draw_ascii(tree)
