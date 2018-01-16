import os
import numpy as np
import igraph as ig

def read_graphs_setup_matrix(root, subset):
    graphs = []
    for patient_id in subset or os.listdir(root):
        graphs.append(ig.Graph.Read_GraphML(os.path.join(root, patient_id, 'optimal.graphml')))
    n = len(graphs)
    matrix = np.zeros((n, n))
    return graphs, matrix

def populate_matrix_symmetric(graphs, matrix, metric, **kwargs):
    for i, G1 in enumerate(graphs):
        for j, G2 in enumerate(graphs[:i+1]):
            val = metric(G1, G2, **kwargs)
            matrix[i][j], matrix[j][i] = val, val
    return matrix

def populate_matrix(graphs, matrix, metric, **kwargs):
    for i, G1 in enumerate(graphs):
        for j, G2 in enumerate(graphs):
            val = metric(G1, G2, **kwargs)
            matrix[i][j], matrix[j][i] = val, val
    return matrix

def node_intersection_count(G1, G2):
    '''

    '''
    nodes1 = { v['name'] for v in G1.vs }
    return len(nodes1.intersection({v['name'] for v in G2.vs}))

def node_intersection_count_matrix(root, subset=None):
    '''

    '''
    graphs, matrix = read_graphs_setup_matrix(root, subset)
    return populate_matrix_symmetric(graphs, matrix, node_intersection_count)

# TODO: exclude node overlap

def shortest_path_value_average(G1, G2, shortest_paths):
    return sum(shortest_paths[(u['name'],v['name'])] for u in G1.vs for v in G2.vs) / (len(G1.vs) * len(G2.vs)) / 2

def shortest_path_value_max(G1, G2, shortest_paths):
    return max(shortest_paths[(u['name'], v['name'])] for u in G1.vs for v in G2.vs)

def shortest_path_value_min(G1, G2, shortest_paths):
    return min(shortest_paths[(u['name'], v['name'])] for u in G1.vs for v in G2.vs)

def _shortest_paths(nodes, graph, **kwargs):
    vs = graph.vs.select(name_in=nodes)
    sps = graph.shortest_paths_dijkstra(source=vs, target=vs, **kwargs)
    return { (u['name'],v['name']): sps[i][j] for i, u in enumerate(vs) for j, v in enumerate(vs) }

def shortest_path_value_matrix(root, graph, subset=None, strategy='average', **kwargs):
    graphs, matrix = read_graphs_setup_matrix(root, subset)
    nodes = {v['name'] for G in graphs for v in G.vs}
    sps = _shortest_paths(nodes, graph, **kwargs)
    return populate_matrix(graphs, matrix, eval('shortest_path_value_'+strategy), shortest_paths=sps)

