'''

'''
import os

import requests
import pandas as pd
import igraph as ig

from biomap import BioMap

################################################################################
# data path                                                                    #
################################################################################

DEREGNET_GRAPH_DATA = os.path.expanduser('~/.deregnet/graphs')
if not os.path.isdir(DEREGNET_GRAPH_DATA):
    os.makedirs(DEREGNET_GRAPH_DATA)

################################################################################
# functionality for parsing edge tables to graphs                              #
################################################################################

def table_to_graph(table,
                   source_column,
                   target_column,
                   directed=True,
                   attributes=None,
                   make_edges_unique=True,
                   exclude=None,
                   include=None,
                   graph_attributes=None,
                   **kwargs):

    return ig.Graph(**table_to_igraph_init_kwargs(table,
                                                  source_column,
                                                  target_column,
                                                  directed,
                                                  attributes,
                                                  make_edges_unique,
                                                  exlcude,
                                                  include,
                                                  graph_attributes,
                                                  **kwargs))
    #

def table_to_igraph_init_kwargs(table,
                                source_column,
                                target_column,
                                directed=True,
                                attributes=None,
                                make_edges_unique=True,
                                exclude=None,
                                include=None,
                                graph_attributes=None,
                                **kwargs):

    def handle_entrez_like_ids(ID):
        try:
            ID = str(int(float(ID)))
        except:
            ID = ID
        return ID

    #
    attributes = {} if attributes is None else attributes
    graph_attributes = {} if graph_attributes is None else graph_attributes
    #
    #
    if isinstance(table, pd.DataFrame):
        data = table
    else:
        data = pd.read_table(table, **kwargs)
    # drop NA's
    data.dropna(inplace=True)
    # filter
    if exclude:
        for fltr in exclude:
            data = data.loc[~(data[fltr['attr']].isin(fltr['values']))]
    elif include:
        for fltr in include:
            data = data.loc[data[fltr['attr']].isin(fltr['values'])]
    # 
    # return value
    kwargs = {}
    kwargs['directed'] = directed
    kwargs['graph_attrs'] = graph_attributes
    # nodes
    nodes = set(data[source_column]) | set(data[target_column])
    nodes = [handle_entrez_like_ids(node) for node in nodes if node.strip()]
    node_index = {node: nodes.index(node) for node in nodes}
    kwargs['n'] = len(nodes)
    kwargs['vertex_attrs'] = {'name': nodes}
    # edges
    edge_attrs = {}
    data['edges'] = list(zip(data[source_column], data[target_column]))
    if make_edges_unique:
        for attr in attributes:
            edge_attrs[attr] = data.groupby('edges')[attr].apply(list).to_dict()
        data.drop_duplicates('edges', inplace=True)
    else:
        for attr in attributes:
            edge_attrs[attr] = data[attr].tolist()
    edges = data['edges'].tolist()
    edge_attrs = {
                   attributes[attr]: [list(set(edge_attrs[attr][edge])) for edge in edges]
                                     if make_edges_unique else edge_attrs[attr]
                   for attr in attributes
                 }
    edges = [
              (node_index[edge[0]],
               node_index[edge[1]])

              for edge in edges
            ]
    kwargs['edges'] = edges
    kwargs['edge_attrs'] = edge_attrs
    return kwargs

def read_sif(sif,
             directed=True,
             make_edges_unique=True,
             exclude=None,
             include=None,
             **kwargs):
    '''
    SIF (Simple Interaction Format) Reader.

    http://wiki.cytoscape.org/Cytoscape_User_Manual/Network_Formats

    WARNING: This function does not support the "a pp x,y,z" notation for
    multiple edges in one line. Instead, a line like the one above will
    be interpreted as one edge between a node "a" and a node "x,y,z".

    Args:

        path (str) : Path of sif file
        directed (bool) : Whether to interpret the graph as directed

    Returns:

        ig.Graph : Graph encoded in the SIF file (hopefully;)
    '''
    return ig.Graph(**sif_to_igraph_init_kwargs(sif,
                                                directed,
                                                make_edges_unique,
                                                exclude,
                                                include,
                                                **kwargs))

def sif_to_igraph_init_kwargs(sif,
                              directed=True,
                              make_edges_unique=True,
                              exclude=None,
                              include=None,
                              **kwargs):
    '''

    '''
    exclude = [{'attr': 1, 'values': exclude}] if exclude is not None else None
    include = [{'attr': 1, 'values': include}] if include is not None else None
    return table_to_igraph_init_kwargs(table=sif,
                                       source_column=0,
                                       target_column=2,
                                       directed=directed,
                                       attributes={1:'interactions' if make_edges_unique else 'interaction'},
                                       make_edges_unique=make_edges_unique,
                                       exclude=exclude,
                                       include=include,
                                       header=None,
                                       **kwargs)

################################################################################
# DeregnetGraph                                                                #
################################################################################

class DeregnetGraph(ig.Graph):

    def __init__(self,
                 make_edges_unique=True,
                 **igraph_init_kwargs):
        '''

        '''
        self.list_valued_interaction_type_attr = make_edges_unique
        super().__init__(**igraph_init_kwargs)

    @property
    def undirected_edge_types(self):
        raise NotImplementedError

    @property
    def edge_type_attribute(self):
        return 'interactions' if self.list_valued_interaction_type_attr else 'interaction'

    def _download(self, url, local_file, verbose):
        if verbose:
            print('Downloading %s ...' % url)
        response = requests.get(url)
        if response.status_code != 200:
            print('Download of %s FAILED.' % url)
        with open(local_file, 'wb') as fp:
            fp.write(response.content)

    def download(self, *args, **kwargs):
        raise NotImplementedError

    def map_nodes(self, mapper, FROM=None, TO=None, source_attr='name', target_attr='name'):
        self.vs[target_attr] = mapper.map(self.vs[source_attr], FROM, TO)

    def map_nodes_to_multiple_targets(self, mapper, FROM=None, TO=[None], source_attr='name', target_attrs=['name']):
        for target_id_type, target_attr in zip(TO, target_attrs):
            self.map_nodes(mapper, FROM, target_id_type, source_attr, target_attr)

    def map_nodes_from_dict(self, dct, source_attr='name', target_attr='name'):
        pass

    def change_name_attr(self, new_name_attr, old_name_attr='_name'):
        self.vs[old_name_attr] = list(self.vs['name'])
        self.vs['name'] = list(self.vs[new_name_attr])

    @property
    def interaction_types(self):
        if self.list_valued_interaction_type_attr:
            return { interaction_type for edge in self.es
                                      for interaction_type in edge[self.edge_type_attribute] }
        else:
            return { edge[self.edge_type_attribute] for edge in self.es }

    def expand_nodes(self, node_attr, keep):
        '''

        '''
        node_index = 0
        names = []
        oldidx2index = {}
        index2targets = {}
        index2sources = {}
        for node in self.vs:
            if node[node_attr] is None:
                key = keep[0]
                value = keep[1]
                if not node[key] == value:
                    continue
                names.append(node['name'])
                oldidx2index[node.index] = node_index
                node_index += 1
            else:
                for attr_val in node[node_attr]:
                    names.append(attr_val)
                    oldidx2index[node.index] = node_index
                    node_index += 1
        num_nodes = node_index
        edge2attrs = {(oldidx2index.get(edge.source, None), oldidx2index.get(edge.target, None)): edge.attributes() for edge in self.es}
        edge2attrs = { edge: edge2attrs[edge] for edge in edge2attrs if edge[0] is not None and edge[1] is not None }
        edges = list(edge2attrs.keys())
        edge_attrs = {}
        for attr in self.es.attribute_names():
            edge_attrs[attr] = [edge2attrs[edge][attr] for edge in edges]

        return DeregnetGraph(self.list_valued_interaction_type_attr,
                             **{'n': num_nodes,
                                'edges': edges,
                                'vertex_attrs': {'name': names},
                                'edge_attrs': edge_attrs,
                                'directed': self.is_directed()})


    def direct_undirected_edges(self, is_directed):
        pass

################################################################################
# Reactome FI 
################################################################################


class ReactomeFI(DeregnetGraph):

    REACTOME_FI_DOWNLOAD_URL = 'http://reactomews.oicr.on.ca:8080/caBigR3WebApp2016'
    FILENAME = 'FIsInGene_022717_with_annotations.txt.zip'

    def __init__(self, exclude=None, include=None, direct_undirected=False, verbose=True, local_root=None):
        self.verbose = verbose
        self.download_root = self.REACTOME_FI_DOWNLOAD_URL
        self.local_root = os.path.join(DEREGNET_GRAPH_DATA, 'reacomte_fi') if local_root is None else local_root
        if not os.path.isdir(self.local_root):
            os.makedirs(self.local_root)
        if not os.path.isfile(self.filepath):
            self.download()
        # only directed edges
        exclude = [{'attr': 'Direction', 'values': ['-']}] if exclude is None else exclude
        igraph_init_kwargs = table_to_igraph_init_kwargs(table=self.filepath,
                                                         source_column='Gene1',
                                                         target_column='Gene2',
                                                         directed=True,
                                                         attributes={
                                                                      'Annotation': 'interaction',
                                                                      'Direction': 'direction',
                                                                      'Score': 'score'
                                                                    },
                                                         make_edges_unique=False,
                                                         exclude=exclude,
                                                         include=include,
                                                         compression='zip')
        super().__init__(make_edges_unique=False,
                         **igraph_init_kwargs)


        if direct_undirected:
            self.direct_undirected_edges(lambda e: e['direction'] == '-')

        self.es['interactions'] = [
                                    [interaction.strip() for interaction in edge['interaction'].split(';')]
                                    for edge in self.es
                                  ]
        self.list_valued_interaction_type_attr = True

        self.vs['symbol'] = self.vs['name']
        self.map_nodes_to_multiple_targets(BioMap.get_mapper('hgnc'),
                                           TO=['entrez', 'ensembl', 'uniprot_ids', 'mgi_id'],
                                           target_attrs=['entrez', 'ensembl', 'uniprot_ids', 'mgi_id'])

    def map_to_mouse(self):
        self.expand_nodes('mgi_id')

    @property
    def filepath(self):
        return os.path.join(self.local_root, self.FILENAME)

    def download(self):
        url = self.download_root+'/'+self.FILENAME
        self._download(url, self.filepath, self.verbose)

    @property
    def undirected_edge_types(self):
        return {edge['interaction'] for edge in self.es if edge['direction'] == '-'}


################################################################################
# KEGG                                                                         #
################################################################################

class KEGG(DeregnetGraph):

    def __init__(self, species='hsa', make_edges_unique=True):
        super().__init__(make_edges_unique, **self._igraph_init(species))

    @classmethod
    @property
    def undirected_edge_types(cls):
        return {'compound',
                'binding/association',
                'dissociation',
                'missing interaction',
                'NA'}

    @property
    def edge_type_attribute(self):
        return 'interactions' if self.make_edges_unique else 'interaction'

    @classmethod
    def download(self, species='hsa'):
        pass

    @classmethod
    def get(self, species='hsa'):
        pass

    @classmethod
    def _igraph_init(self, species='hsa'):
        pass

    # TODO: implement download


################################################################################
# Omnipath 
################################################################################

class OmniPath(DeregnetGraph):
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

class PathwayCommons(DeregnetGraph):

    def __init__(self,
                 what='All',
                 exclude=[],
                 include=[],
                 make_edges_unique=True,
                 download_root=PATHWAY_COMMONS_DOWNLOAD_ROOT,
                 version=9,
                 verbose=True):

        self.root = download_root
        self.version = version
        self.verbose = verbose
        if not os.path.isdir(self.local_path):
            os.makedirs(self.local_path)
        filename = self._file_name(what)
        filepath = os.path.join(self.local_path, filename)
        if not os.path.isfile(filepath):
            self.download(what)

        igraph_init_kwargs = sif_to_igraph_init_kwargs(sif=filepath,
                                                       directed=True,
                                                       make_edges_unique=make_edges_unique,
                                                       exclude=exclude,
                                                       include=include,
                                                       compression='gzip')

        super().__init__(make_edges_unique, **igraph_init_kwargs)

        hgnc = BioMap.get_mapper('hgnc')
        self.vs['symbol'] = hgnc.map(self.vs['name'], TO='symbol')
        self.map_nodes_to_multiple_targets(hgnc,
                                           TO=['entrez', 'ensembl', 'uniprot_ids', 'mgd_id'],
                                           target_attrs=['entrez', 'ensembl', 'uniprot_ids', 'mgi_id'])

    def map_to_mouse(self):
        mouse_graph = self.expand_nodes('mgi_id', keep=('symbol', None))
        mgi_ensembl = BioMap.get_mapper('mgi_ensembl')
        mouse_graph.map_nodes_to_multiple_targets(mgi_ensembl, TO=['symbol', 'ensembl_id'], target_attrs=['symbol', 'ensembl'])
        mgi_entrez = BioMap.get_mapper('mgi_entrez')
        mouse_graph.map_nodes(mgi_entrez, TO='entrez_id', target_attr='entrez')
        return mouse_graph

    def download(self, what):
        filename = self._file_name(what)
        url = self._download_url(filename)
        filepath = os.path.join(self.local_path, filename)
        self._download(url, filepath, self.verbose)

    @property
    def local_path(self):
        return os.path.join(DEREGNET_GRAPH_DATA, 'pathway_commons')

    def _download_url(self, filename):
        return self.root+'/v'+str(self.version)+'/'+filename

    def _file_name(self, what):
        return 'PathwayCommons'+str(self.version)+'.'+what+'.hgnc.sif.gz'

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


    def download_all(self):
        for data_source in self.available_data_sources:
            self.download(data_source)

    @classmethod
    def undirected_edge_types(cls):
        return { 'reacts-with',
                 'interacts-with',
                 'in-complex-with' }

