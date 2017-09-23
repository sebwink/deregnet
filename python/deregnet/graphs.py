import os
import re
import gzip
import zipfile

import requests
import igraph as ig

DEREGNET_GRAPH_DATA=os.path.expanduser('~/.deregnet/graphs')

if not os.path.isdir(DEREGNET_GRAPH_DATA):
    os.makedirs(DEREGNET_GRAPH_DATA)

def tostr(bytes_or_str):
    if isinstance(bytes_or_str, bytes):
        return bytes_or_str.decode('utf-8')
    return bytes_or_str

def read_sif(sif, directed=True):
    '''
    SIF (Simple Interaction Format) Reader.

    http://wiki.cytoscape.org/Cytoscape_User_Manual/Network_Formats

    Args:

        path (str) : Path of sif file
        directed (bool) : Whether to interpret the graph as directed

    Returns:

        ig.Graph : Graph encoded in the SIF file (hopefully;)
    '''
    nodes = set()
    edges = list()
    interactions = list()
    for line in sif:
        line = tostr(line)
        items = [item.strip() for item in re.split('\s+', line) if item]
        if len(items) == 1:    # isolated nodes
            nodes.add(source)
            continue
        source, interaction, targets = items[0], items[1], items[2:]
        nodes.add(source)
        nodes |= {target for target in targets}
        edges.extend([(source, target) for target in targets])
        interactions.extend( len(targets) * [interaction] )
    nodes = list(nodes)
    node2index = { node : nodes.index(node) for node in nodes }
    edges = [(node2index[edge[0]], node2index[edge[1]]) for edge in edges]
    graph = ig.Graph(directed=directed)
    graph.add_vertices(len(nodes))
    graph.vs['name'] = nodes
    graph.add_edges(edges)
    graph.es['interaction'] = interactions
    return graph

def map_node_ids(graph, id_mapper, id_attr, id_type, id_targets):
    pass

################################################################################
# KEGG 
################################################################################

KEGG_UNDIRECTED_EDGE_TYPES = { 'compound',
                               'binding/association',
                               'dissociation',
                               'missing interaction',
                               'NA' }

class KEGG:

    # TODO: implement download

    undirected_edge_types = KEGG_UNDIRECTED_EDGE_TYPES

    def __init__(self, path_or_species='hsa'):
        if path_or_species in {'hsa', 'mmu', 'rno', 'sce'}:
            species = path_or_species
            self.path = os.path.join(DEREGNET_GRAPH_DATA, 'kegg/kegg_'+species+'.graphml')
        else:
            self.path = path_or_species

    def __call__(self, exclude=KEGG_UNDIRECTED_EDGE_TYPES, include=None):
        graph = ig.Graph.Read_GraphML(self.path)
        graph.es['interaction'] = [interaction.split(',') for interaction in graph.es['interaction']]
        if exclude:
            self.exclude_edge_types(graph, exclude)
        if (not exclude) and include:
            self.include_edge_types(graph, include)
        # direct undirected edge types
        for edge in graph.es:
            for interaction in self.undirected_edge_types:
                self._direct_undirected_edge_types(graph, interaction, edge, edge['interaction'])
        return graph

    def exclude_edge_types(self, graph, interactions):
        interactions = set(interactions)
        edges_to_delete = set()
        for edge in graph.es:
            diff = set(edge['interaction']) - interactions
            if not diff:
                edges_to_delete.add(edge.index)
            else:
                edge['interaction'] = list(diff)
        graph.delete_edges(edges_to_delete)

    def include_edge_types(self, graph, interactions):
        interactions = set(interactions)
        edges_to_delete = set()
        for edge in graph.es:
            intersection = set(edge['interaction']) & interactions
            if not intersection:
                edges_to_delete.add(edge.index)
            else:
                edge['interaction'] = list(intersection)
        graph.delete_edges(edges_to_delete)

    def _direct_undirected_edge_types(self, graph, interaction, edge, edge_interactions):
        if interaction in edge_interactions:
            rev_edge = graph.es.select(_between=({edge.target},{edge.source}))
            if len(rev_edge) > 0:
                list(rev_edge)[0]['interaction'] = list(set(edge_interactions) | {interaction})
            else:
                graph.add_edge(edge.target,
                               edge.source,
                               interaction=interaction)

################################################################################
# Omnipath 
################################################################################

class OmniPath:
    # TODO: implement download here

    def __init__(self, path=None):
        if path is None:
            self.path = DEREGNET_GRAPH_DATA
        else:
            self.path = path

    def __call__(self):
        return ig.Graph.Read_GraphML(os.path.join(self.path, 'omnipath/omnipath_directed_interactions.graphml'))

    def ptm_graph(self):
        return ig.Graph.Read_GraphML(os.path.join(self.path, 'omnipath/omnipath_ptm_graph.graphml'))


################################################################################
# Pathway Commons
################################################################################

PATHWAY_COMMONS_DOWNLOAD_ROOT='http://www.pathwaycommons.org/archives/PC2'

class PathwayCommons:

    def __init__(self, download_root=PATHWAY_COMMONS_DOWNLOAD_ROOT, version=9, verbose=True):
        self.root = download_root
        self.local_path = os.path.join(DEREGNET_GRAPH_DATA, 'pathway_commons')
        if not os.path.isdir(self.local_path):
            os.makedirs(self.local_path)
        self.version = version
        self.verbose = verbose

    def _download(self, filename):
        url = self.root+'/v'+str(self.version)+'/'+filename
        if self.verbose:
            print('Downloading %s ...' % url)
        response = requests.get(url)
        if response.status_code != 200:
            print('Download of %s failed.' % url)
            return None
        local_file = os.path.join(self.local_path, filename)
        with open(local_file, 'wb') as fp:
            fp.write(response.content)

    @property
    def available_data_sources(self):
        return {
                 'wp': '',
                 'smpdb': '',
                 'reconx': '',
                 'reactome': '',
                 'psp': '',
                 'pid': '',
                 'panther': '',
                 'netpath': '',
                 'msigdb': '',
                 'kegg': '',
                 'intact': '',
                 'intact_complex': '',
                 'inoh': '',
                 'humancyc': '',
                 'hprd': '',
                 'drugbank': '',
                 'dip': '',
                 'ctd': '',
                 'corum': '',
                 'biogrid': '',
                 'bind': '',
                 'All': '',
                 'Detailed': ''
               }


    def download(self, what):
        filename = 'PathwayCommons'+str(self.version)+'.'+what+'.hgnc.sif.gz'
        self._download(filename)

    def download_all(self):
        for data_source in self.available_data_sources:
            self.download(data_source)

    def _get(self, path):
        with gzip.open(path, 'rb') as sif:
            return read_sif(sif)

    def get(self, what):
        filename = 'PathwayCommons'+str(self.version)+'.'+what+'.hgnc.sif.gz'
        filepath = os.path.join(self.local_path, filename)
        if not os.path.isfile(filepath):
            self.download(what)
        graph = self._get(filepath)
        if self.verbose:
            print('Requested base graph \'%s\': %s nodes, %s edges' % (what, str(len(graph.vs)), str(len(graph.es))))
        return graph

    def __call__(self, what, exclude=None, include=None):
        graph = self.get(what)
        if exclude:
            self.exclude_interaction_types(graph, exlcude)
        if not exclude and include:
            self.include_interactions_types(graph, include)
        if self.verbose:
            print('Returning filtered graph: %s nodes, %s edges' % (str(len(graph.vs)), str(len(graph.es))))
        return graph

    def exclude_interaction_types(self, graph, exclude):
        edges_to_delete = {edge for edge in graph.es if edge['interaction'] in exclude}
        graph.delete_edges(edges_to_delete)

    def include_interactions_types(self, graph, include):
        edges_to_delete = {edge for edge in graph.es if edge['interaction'] not in include}
        graph.delete_edges(edges_to_delete)

################################################################################
# Reactome FI 
################################################################################

REACTOME_FI_DOWNLOAD_URL = 'http://reactomews.oicr.on.ca:8080/caBigR3WebApp2016'

class ReactomeFI:

    FILENAME = 'FIsInGene_022717_with_annotations.txt.zip'

    def __init__(self, download_url=REACTOME_FI_DOWNLOAD_URL):
        self.root = download_url
        self.local_path = os.path.join(DEREGNET_GRAPH_DATA, 'reactome_fi')
        if not os.path.isdir(self.local_path):
            os.makedirs(self.local_path)

    def _download(self, filename):
        url = self.root+'/'+filename
        response = requests.get(url)
        if response.status_code != 200:
            print('Download FAILED.')
        local_file = os.path.join(self.local_path, filename)
        with open(local_file, 'wb') as fp:
            fp.write(response.content)

    def download(self):
        self._download(self.FILENAME)

    def parse(self):
        local_file = os.path.join(self.local_path, self.FILENAME)
        with zipfile.ZipFile(local_file) as z:
            with z.open('.'.join(self.FILENAME.split('.')[:-1])) as f:
                return [tostr(line)[:-1].split('\t') for line in f]

    def get(self):
        edges = self.parse()[1:]
        for edge in edges:
            pass

