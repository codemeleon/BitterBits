import networkx as nx
from ete3 import Tree
from pylab import *
# %%
tree = Tree(
    "/home/devil/Documents/Influenza/NewAnalysis/Trees/RAxML_bestTree.FluB_iva_NA__aln_trim_select"
)
tree.ladderize(direction=0)
x_dist = {}
for i, node in enumerate(tree.traverse("postorder")):
    if node.name == "NoName" or not node.is_leaf():
        node.name = str(i)
    x_dist[node.name] = -1 * tree.get_distance(node)

# x_core
y_dist = {}
nodes = []
leafs = tree.get_leaf_names()
number_of_leafs = len(leafs)
# nodes = list(tree.get_leaves())
# nodes.reverse()
for j in range(number_of_leafs):
    y_dist[leafs[j]] = number_of_leafs - j
    nodes.append(leafs[j])


def ycoors(tree):  # Calculate Y_distances of nodes from the bottom
    for child in tree.get_children():
        if child.name not in nodes:
            ycoors(child)
    y_dist[tree.name] = (y_dist[tree.get_children()[0].name] +
                         y_dist[tree.get_children()[1].name]) / 2.0
    nodes.append(tree.name)


labels = {}
for k in leafs:
    labels[k] = k

ycoors(tree)

edges = []
int_nodes = [] # configure intemediate node to make tree look good
for node in tree.traverse("preorder"):
    if node.is_leaf():
        continue
    else:
        add_nodes = 1
        if add_nodes:
            c0 = node.get_children()[0].name
            c1 = node.get_children()[1].name
            p = node.name
            p_c0 = "%s_%s" % (p, c0)
            p_c1 = "%s_%s" % (p, c1)
            edges += [
                [p, p_c0],
                [p, p_c1],
                [p_c0, c0],
                [p_c1, c1]
            ]
            x_dist[p_c0] = x_dist[p]
            x_dist[p_c1] = x_dist[p]
            y_dist[p_c0] = y_dist[c0]
            y_dist[p_c1] = y_dist[c1]





            pass
        else:
            edges += [[node.name, node.get_children()[0].name], # need to add the positions of the new nodes
                    [node.name, node.get_children()[1].name]]
nodes_coor = {}
for k in y_dist:
    nodes_coor[k] = (x_dist[k], y_dist[k])

   # print(node.name)
graph = nx.Graph()
graph.add_edges_from(edges)

nx.draw(graph, nodes_coor, node_size=1, labels=labels, fontsize=0.1)
show()
