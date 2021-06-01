# Modules
import click
from ete3 import Tree
import networkx as nx
from pylab import *


def coors_label(tree, reverse=False, gap=0):
    tree = Tree(tree)
    tree.ladderize(direction=0)
    xcoor = {}
    labels = {}
    for i, node in enumerate(tree.traverse("preorder")):
        if node.is_leaf():
            name = node.name
            if reverse:
                node.name += "_1"
            else:
                node.name += "_0"
            labels[node.name] = name
        else:
            if reverse:
                node.name = "1_%d"%i
            else:
                node.name = "0_%d" %i
        xcoor[node.name] = tree.get_distance(node)

    ycoor = {}
    nodes = []
    leafs = tree.get_leaf_names()
    number_of_leafs = len(leafs)
    for j in range(number_of_leafs):
        ycoor[leafs[j]] = number_of_leafs - j
        nodes.append(leafs[j])

    def ycoorfn(tree):
        for child in tree.get_children():
            if child.name not in nodes:
                ycoorfn(child)
        if len(tree.get_children()) > 2:
            print(tree.name)
            for child in tree.get_children():
                print(child.name, tree.get_distance(child))
        ycoor[tree.name] = (ycoor[tree.get_children()[0].name] +
                            ycoor[tree.get_children()[1].name]) / 2.0

        nodes.append(tree.name)
    ycoorfn(tree)
    edges = []
    add_internal_nodes = 1
    for node in tree.traverse("preorder"):
        if node.is_leaf():
            continue
        else:
            if add_internal_nodes:
                if not node.is_root():
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
                    xcoor[p_c0] = xcoor[p]
                    xcoor[p_c1] = xcoor[p]
                    ycoor[p_c0] = ycoor[c0]
                    ycoor[p_c1] = ycoor[c1]
                else:
                    # Note: For root node there are three children. Please
                    # correct it for other situations later
                    # c0 = node.get_children()[0].name
                    c0 = node.get_children()[1].name
                    c1 = node.get_children()[2].name
                    p = node.name
                    p_c0 = "%s_%s" % (p, c0)
                    p_c1 = "%s_%s" % (p, c1)
                    edges += [
                        [p, p_c0],
                        [p, p_c1],
                        [p_c0, c0],
                        [p_c1, c1]
                    ]
                    xcoor[p_c0] = xcoor[p]
                    xcoor[p_c1] = xcoor[p]
                    ycoor[p_c0] = ycoor[c0]
                    ycoor[p_c1] = ycoor[c1]
                # need to add the positions of the new nodes
            else:
                edges += [[node.name, node.get_children()[0].name],
                          [node.name, node.get_children()[1].name]]
    nodes_coor = {}
    if reverse:
        shift = max(xcoor.values()) + gap
        for k in xcoor:
            xcoor[k] = -1*xcoor[k] + shift
    for k in ycoor:
        nodes_coor[k] = (xcoor[k], ycoor[k])
    return edges, nodes_coor, labels, max(xcoor.values())


@click.command()
@click.option(
    "-t1",
    help="First Tree",
    type=str,
    default="/home/devil/Documents/Influenza/NewAnalysis/Trees/RAxML_bestTree.FluB_iva_HA_aln_trim_select",
    required=True,
    show_default=True)
@click.option(
    "-t2",
    help="First Tree",
    type=str,
    default="/home/devil/Documents/Influenza/NewAnalysis/Trees/RAxML_bestTree.FluB_iva_NA__aln_trim_select",
    required=True,
    show_default=True)

def run(t1,t2):
    edges1, nodes_coors1, labels1, max_x = coors_label(t1)
    print(max_x, "Anmol")
    edges2, nodes_coors2, labels2, _ = coors_label(t2, reverse=True,
                                                   gap=max_x)
    nodes_coors1 = {**nodes_coors1, **nodes_coors2}
    labels = {**labels1, **labels2}
    graph = nx.Graph()
    # print(nodes_coors2)
    #exit(1)

    graph.add_edges_from(edges1+edges2)
    # graph.add_edges_from(edges2)
    ## Arrange the coordinates here

    # nx.draw(graph, nodes_coors1, node_size=1,  labels=labels1, fontsize=0.1)
    nx.draw(graph, nodes_coors1, node_size=0)
    # Adding text
    for lb in labels1:
        x,y = nodes_coors1[lb]
        text(x, y, labels1[lb], fontsize=5, ha='left', va='center')
    for lb in labels2:
        x,y = nodes_coors1[lb]
        text(x, y, labels2[lb], fontsize=5, ha='right', va='center')

    # Connected samples
    reverse_labels1 = {}
    for k in labels1:
        reverse_labels1[labels1[k]] = k

    reverse_labels2 = {}
    for k in labels2:
        reverse_labels2[labels2[k]] = k
    common_labels = set(reverse_labels1)&set(reverse_labels2)
    for cl in common_labels:
        x = (nodes_coors1[reverse_labels1[cl]][0],nodes_coors2[reverse_labels2[cl]][0])
        y = (nodes_coors1[reverse_labels1[cl]][1],nodes_coors2[reverse_labels2[cl]][1])
        print(cl,x,y)
        plot(x,y,c='k')
    show()


if __name__ == '__main__':
    run()
