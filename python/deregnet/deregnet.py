import os
import time
import subprocess
import igraph
import numpy as np
import networkx as nx
import pandas as pd


__tmpdir__ = os.path.join(os.path.expanduser('~'), '.deregnet/tmp')

if not os.path.isdir(__tmpdir__):
    os.makedirs(__tmpdir__)

def time_stamp():
    return time.strftime('%Y-%m-%d-%H-%M-%S"', time.gmtime())

def igraph2lgf(graph, path,
                      keep_node_attr = [],
                      keep_edge_attr = [],
                      kepp_graph_attr = [],
                      additional_graph_attributes = {},
                      id_attr = 'id'):
    '''
    Writes an igraph graph to lemon graph format (lgf).
    '''
    def write_node(node, node_label, lgf):
        lgf.write(str(node_label) + '\t' + node[id_attr] + '\n')

    def write_nodes(lgf):
        lgf.write('@nodes\n')
        lgf.write('label\tid\n')
        for node_label, node in enumerate(graph.vs):
            write_node(node, node_label, lgf)

    def write_edge(edge, edge_label, lgf):
        lgf.write(str(edge.source) + '\t' + \
                  str(edge.target) + '\t' + \
                  str(edge_label) + '\n')

    def write_edges(lgf):
        lgf.write('@arcs\n')
        lgf.write('\t\tlabel\n')
        for edge_label, edge in enumerate(graph.es):
            write_edge(edge, edge_label, lgf)

    def write_graph_attributes(attributes, lgf):
        lgf.write('@attributes\n')
        for attr in attributes:
            lgf.write(attr + '\t' + str(attributes[attr]) + '\n')

    with open(path, 'w') as lgf:
        write_nodes(lgf)
        write_edges(lgf)
        if additional_graph_attributes:
            write_graph_attributes(additional_graph_attributes, lgf)

class TmpFileHandle: 
    def __init__(self, tmpdir = __tmpdir__):
        self.path = os.path.join(tmpdir, time_stamp()[:-1])
        os.mkdir(self.path)
        os.mkdir(os.path.join(self.path, 'subgraphs'))
        self.graph = None
        self.scores = None
    def write_graph(self, graph, graph_id_attr):
        self.graph = os.path.join(self.path, 'graph.lgf')
        igraph2lgf(graph, self.graph, id_attr = graph_id_attr)
    def write_scores(self, scores):
        self.scores = os.path.join(self.path, 'scores.tsv')
        scores.to_csv(self.scores, sep = '\t', index = False, header = False)

def parse_scores(path2score,
                 id_col,
                 score_col,
                 sep,
                 score_id_type,
                 graph_id_type,
                 species,
                 id_mapper,
                 **kwargs):
    df = pd.read_csv(path2score, sep = sep)
    df = df[[id_col, score_col]]
    df.dropna(inplace=True)
    ids = list(df[id_col])
    if id_mapper is None:
        graph_id_type = score_id_type
    if score_id_type != graph_id_type:
        mapper = id_mapper(default_species = species, **kwargs)
        ids = mapper.map(ids, score_id_type, graph_id_type)
    df['ids'] = ids
    df = df[['ids', score_col]]
    df = df[df['ids'] != '']
    # df = df[str(df[score_col]) != '']
    return df

def remove_self_loops(graph):
    self_loops = [e for e in graph.es if e.source == e.target]
    graph.delete_edges(self_loops)
    return graph

def write_tmp_files(graph,
                    graph_id_attr,
                    graph_id_type,
                    path2score,
                    score_id_type,
                    sep,
                    id_col,
                    score_col,
                    species,
                    id_mapper = None,
                    **kwargs):
    graph = remove_self_loops(graph)
    scores = parse_scores(path2score,
                          id_col,
                          score_col,
                          sep,
                          score_id_type,
                          graph_id_type,
                          species,
                          id_mapper,
                          **kwargs)
    tmp_files = TmpFileHandle()
    tmp_files.write_graph(graph, graph_id_attr)
    tmp_files.write_scores(scores)
    return tmp_files

# work around igraph graphml bug
def ig_to_nx(ig_graph):
    def series_to_str(series):
        series = [str(elmt) for elmt in series]
        return '['+','.join(series)+']'
    nx_graph = nx.DiGraph()
    for i, node in enumerate(ig_graph.vs):
        nx_graph.add_node(i, node.attributes())
    for edge in ig_graph.es:
        nx_graph.add_edge(edge.source, edge.target, edge.attributes())
    for node in nx_graph.nodes(data=True):
        for attr in node[1]:
            if isinstance(node[1][attr], np.float64):
                node[1][attr] = float(node[1][attr])
            if isinstance(node[1][attr], pd.core.series.Series):
                node[1][attr] = series_to_str(node[1][attr])
    for edge in nx_graph.edges(data=True):
        for attr in edge[2]:
            if isinstance(edge[2][attr], pd.core.series.Series):
                edge[2][attr] = series_to_str(edge[2][attr])
    return nx_graph


def output2graphml(graph, tmp_files, outdir):
    def parse_sif(path2sif):
        nodes = set()
        with open(path2sif, 'r') as sif:
            for edge in sif:
                edge = edge.split('\t')
                nodes.add(edge[0])
                nodes.add(edge[2][:-1])
        return nodes 

    def get_root(path):
        with open(path, 'r') as sif:
            first_line = sif.readline()
            return first_line.split(':')[-1].strip()

    
    scores = pd.read_csv(tmp_files.scores, sep = '\t', header = None)
    scores.set_index(0, inplace = True)
    output = os.path.join(tmp_files.path, 'subgraphs/plain')
    for sif in os.listdir(output):
        root = get_root(os.path.join(output, '..', sif))
        nodes = parse_sif(os.path.join(output, sif))
        mapped_nodes = graph.vs.select(name_in=nodes)
        subgraph = graph.subgraph(mapped_nodes, "create_from_scratch")
        for node in subgraph.vs:
            try:
                node['score'] = scores.ix[node['name'],1]
            except:
                node['score'] = 0.0
            if node['name'] == root:
                node['root'] = True
            else:
                node['root'] = False
        graphml = os.path.join(outdir, sif.split('.')[0] + '.graphml')
        subgraph = ig_to_nx(subgraph)
        nx.write_graphml(subgraph, graphml)

