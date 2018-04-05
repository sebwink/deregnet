import os
import numpy as np
import igraph as ig


def read_graphs(root, subset):
    graphs = []
    patients = []
    for patient_id in subset or os.listdir(root):
        if patient_id == 'failed.txt': continue
        subgrph = os.path.join(root, patient_id, 'optimal.graphml')
        if not os.path.isfile(subgrph):
            continue
        graphs.append(ig.Graph.Read_GraphML(subgrph))
        patients.append(patient_id)
    return graphs, patients

def read_graphs_setup_matrix(root, subset):
    graphs, patients = read_graphs(root, subset)
    n = len(graphs)
    matrix = np.zeros((n, n))
    return graphs, matrix, patients

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

def node_intersection_count(G1, G2, scale):
    '''

    '''
    nodes1 = { v['name'] for v in G1.vs }
    nodes2 = { v['name'] for v in G2.vs }
    intersection = len(nodes1.intersection(nodes2))
    if scale:
        return intersection / (len(nodes1) + len(nodes2) - intersection)
    else:
        return intersection

def node_intersection_count_matrix(root, subset=None, scale=True):
    '''

    '''
    graphs, matrix, patients = read_graphs_setup_matrix(root, subset)
    return populate_matrix_symmetric(graphs, matrix, node_intersection_count, scale=scale), patients

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
    graphs, matrix, patients = read_graphs_setup_matrix(root, subset)
    nodes = {v['name'] for G in graphs for v in G.vs}
    sps = _shortest_paths(nodes, graph, **kwargs)
    return populate_matrix(graphs, matrix, eval('shortest_path_value_'+strategy), shortest_paths=sps), patients

#
#

def _boolean_attr_overlap(G1, G2, attr, scale):
    pos1 = {node['name'] for node in G1.vs if node[attr]}
    pos2 = {node['name'] for node in G2.vs if node[attr]}
    nodes1 = {node['name'] for node in G1.vs}
    nodes2 = {node['name'] for node in G2.vs}
    intersection = len(nodes1.intersection(nodes2))
    if scale:
        return len(pos1.intersection(pos2)) / (len(nodes1) + len(nodes2) - intersection)
    else:
        return len(pos1.intersection(pos2)) 

def boolean_attr_overlap_count(root, attr, subset=None, scale=True):
    graphs, matrix, patients = read_graphs_setup_matrix(root, subset)
    return populate_matrix_symmetric(graphs, matrix, _boolean_attr_overlap, attr=attr, scale=scale), patients

def receptor_intersection_count_matrix(root, subset=None, scale=True):
    return boolean_attr_overlap_count(root, attr='deregnet_receptor', subset=subset, scale=scale)

def terminal_intersection_count_matrix(root, subset=None, scale=True):
    return boolean_attr_overlap_count(root, attr='deregnet_terminal', subset=subset, scale=scale)


#
#

def get_vogelstein(path='/home/sebastian/wrk/deregnet/deregnet/data/cancer_genes_vogelstein.txt'):
    genes = set()
    with open(path, 'r') as fp:
        for line in fp.readlines():
            genes.add(line.strip())
    return genes

def genes_in(G, geneset):
    return { node['symbol'] for node in G.vs if node['symbol'] in geneset }

def subgraph2vogelstein(root, subset=None):
    graphs, patients = read_graphs(root, subset)
    vogelstein = get_vogelstein()
    return { patients[i]: genes_in(G, vogelstein) for i, G in enumerate(graphs) }
