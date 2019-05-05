import os
import sys
import igraph as ig

def _consensus_subgraph(graphs, graph):
    nodes = {v['name']: 0 for G in graphs for v in G.vs}
    for G in graphs:
        for v in G.vs:
            nodes[v['name']] += 1
    subgraph = graph.subgraph(graph.vs.select(name_in=nodes.keys()))
    for v in subgraph.vs:
        v['occurence_count'] = nodes[v['name']]
    return subgraph


def consensus_subgraph(root, graph, subset=None):
    graphs = []
    for patient_id in subset or os.listdir(root):
        graphs.append(ig.Graph.Read_GraphML(os.path.join(root, patient_id, 'optimal.graphml')))
    return _consensus_subgraph(graphs, graph)

if __name__ == '__main__':
    graph = ig.Graph.Read_GraphML(sys.argv[2])
    subgraph = consensus_subgraph(sys.argv[1], graph)
    subgraph.write_graphml('consensus_subgraph.graphml')
