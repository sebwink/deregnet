import os
import sys

import igraph as ig

def read_sif(path, directed=True):
    '''

    '''
    nodes = set()
    edges = set()
    interactions = dict()
    with open(path, 'r') as sif:
        for line in sif:
            items = [item.strip() for item in line.split('\t') if item]
            print(items)
            source, interaction, target = items[0].split(':')[1], items[1], items[2].split(':')[1]
            nodes.add(source)
            nodes.add(target)
            edge = (source, target)
            if edge not in interactions:
                interactions[edge] = {interaction}
            else:
                interactions[edge].add(interaction)
    for edge in interactions:
        interactions[edge] = ','.join(list(interactions[edge]))
    nodes = list(nodes)
    print(nodes)
    node2index = { node : nodes.index(node) for node in nodes }
    print(node2index)
    edges = [(node2index[edge[0]], node2index[edge[1]]) for edge in interactions]
    print(edges)
    interactions = [interactions[(nodes[edge[0]], nodes[edge[1]])] for edge in edges]
    print(interactions)
    graph = ig.Graph(directed=directed)
    graph.add_vertices(len(nodes))
    graph.vs['name'] = nodes
    graph.add_edges(edges)
    graph.es['interaction'] = interactions
    return graph

def to_graphml(sif):
    graph = read_sif(sif, True)
    graphml = os.path.join(os.path.dirname(sif), 
                           os.path.basename(sif).split('.')[0]+'.graphml')
    graph.write_graphml(graphml)

if __name__ == '__main__':
    to_graphml(sys.argv[1])
