'''

'''
import os
import zipfile

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
    nodes = [handle_entrez_like_ids(node) for node in nodes]
    nodes = [node for node in nodes if node.strip()]
    node_index = {node: nodes.index(node) for node in nodes}
    kwargs['n'] = len(nodes)
    kwargs['vertex_attrs'] = {'name': nodes}
    # edges
    edge_attrs = {}
    data['edges'] = list(zip(data[source_column].tolist(), data[target_column].tolist()))
    data['edges'] = [(handle_entrez_like_ids(edge[0]), handle_entrez_like_ids(edge[1]))
                     for edge in data['edges']]
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
                                                ** kwargs))

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

    @classmethod
    def _download(cls, url, local_file, verbose):
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

    def map_nodes_from_dict(self, dct, source_attr='name', target_attr='name', default=None):
        self.vs[target_attr] = [dct.get(v[source_attr], default) for v in self.vs]

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

    def neighborhood_graph(self, nodes, mode=ig.ALL, depth=1, node_attr='symbol'):
        # TODO filter
        if isinstance(nodes, str):
            nodes = [nodes]
        neighborhood = set(self.vs.select(**{node_attr+'_in':nodes}))
        for d in range(depth):
            for node in list(neighborhood):
                neighbors = set(self.vs.select(self.neighbors(node)))
                neighborhood =  neighborhood | neighbors
        return self.subgraph(neighborhood, 'create_from_scratch')


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
        self.map_nodes_to_multiple_targets(BioMap().get_mapper('hgnc'),
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

    KEGG_GRAPH_PATH = os.path.join(DEREGNET_GRAPH_DATA, 'kegg')

    def __init__(self,
                 species='hsa',
                 exclude=None,
                 include=None,
                 make_edges_unique=True,
                 directed=True):
        if not os.path.isdir(self.KEGG_GRAPH_PATH):
            os.makedirs(self.KEGG_GRAPH_PATH)
        local_file = os.path.join(self.KEGG_GRAPH_PATH, 'kegg_'+species+'.sif')
        # TODO: implement download
        igraph_init_kwargs = sif_to_igraph_init_kwargs(local_file,
                                                       directed=directed,
                                                       make_edges_unique=make_edges_unique,
                                                       exclude = exclude,
                                                       include = include)

        super().__init__(make_edges_unique, **igraph_init_kwargs)
        self.vs['name'] = [ID.split(':')[-1] for ID in self.vs['name']]
        if species == 'hsa':
            self.map_nodes_to_multiple_targets(BioMap().get_mapper('hgnc'),
                                               FROM='entrez',
                                               TO=['entrez',
                                                   'ensembl',
                                                   'symbol',
                                                   'uniprot_ids'],
                                               target_attrs=['entrez',
                                                             'ensembl',
                                                             'symbol',
                                                             'uniprot_ids'])
        elif species == 'mmu':
            self.map_nodes_to_multiple_targets(BioMap().get_mapper('mgi_entrez'),
                                               FROM='entrez',
                                               TO=['entrez', 'symbol', 'name'],
                                               target_attrs=['entrez', 'symbol', 'name'])
            self.map_nodes(BioMap().get_mapper('mgi_ensembl'),
                           FROM='symbol', TO='ensembl',
                           source_attr='symbol', target_attr='ensembl')

    @classmethod
    def undirected_edge_types(cls):
        return {'binding/association',
                'dissociation',
                'missing interaction',
                'NA'}

    @classmethod
    def download(self, species='hsa'):
        pass

    @classmethod
    def get(self, species='hsa'):
        pass

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

        hgnc = BioMap().get_mapper('hgnc')
        self.vs['symbol'] = hgnc.map(self.vs['name'], TO='symbol')
        self.map_nodes_to_multiple_targets(hgnc,
                                           TO=['entrez', 'ensembl', 'uniprot_ids', 'mgd_id'],
                                            target_attrs=['entrez', 'ensembl', 'uniprot_ids', 'mgi_id'])

    def map_to_mouse(self):
        mouse_graph = self.expand_nodes('mgi_id', keep=('symbol', None))
        mgi_ensembl = BioMap().get_mapper('mgi_ensembl')
        mouse_graph.map_nodes_to_multiple_targets(mgi_ensembl, TO=['symbol', 'ensembl_id'], target_attrs=['symbol', 'ensembl'])
        mgi_entrez = BioMap().get_mapper('mgi_entrez')
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

################################################################################
# RegNetwork #
################################################################################

DEFAULT_REG_NETWORK_DOWNLOAD_URL='http://www.regnetworkweb.org/download'

DEFAULT_REG_NETWORK_LOCAL_PATH=os.path.join(DEREGNET_GRAPH_DATA, 'regnetwork')
if not os.path.isdir(DEFAULT_REG_NETWORK_LOCAL_PATH):
    os.makedirs(DEFAULT_REG_NETWORK_LOCAL_PATH)

class RegNetwork(DeregnetGraph):
    '''
    Gene-regulatory graphs defined for Human and Mouse available at:
     _______________________________________
    |                                       |
    | http://www.regnetworkweb.org/home.jsp |
    |_______________________________________|

    If you use these graphs in your work, please cite (see also RegNetwork.cite):

    -------------------------------------------------------------------------------------------
    Liu et al.: RegNetwork: an integrated database of transcriptional and post-transcriptional
    regulatory networks in human and mouse. Database, 2015, 1-12, doi:10.1093/database/bav095
    -------------------------------------------------------------------------------------------

    Abstract:

    Transcriptional and post-transcriptional regulation of gene expression is of fundamental im-
    portance to numerous biological processes. Nowadays, an increasing amount of gene regu-
    latory relationships have been documented in various databases and literature. However, to
    more efficiently exploit such knowledge for biomedical research and applications, it is ne-
    cessary to construct a genome-wide regulatory network database to integrate the informa-
    tion on gene regulatory relationships that are widely scattered in many different places.
    Therefore, in this work, we build a knowledge-based database, named ‘RegNetwork’, of
    gene regulatory networks for human and mouse by collecting and integrating the docu-
    mented regulatory interactions among transcription factors (TFs), microRNAs (miRNAs) and
    target genes from 25 selected databases. Moreover, we also inferred and incorporated po-
    tential regulatory relationships based on transcription factor binding site (TFBS) motifs into
    RegNetwork. As a result, RegNetwork contains a comprehensive set of experimentally
    observed or predicted transcriptional and post-transcriptional regulatory relationships, and
    the database framework is flexibly designed for potential extensions to include gene regula-
    tory networks for other organisms in the future. Based on RegNetwork, we characterized the
    statistical and topological properties of genome-wide regulatory networks for human and
    mouse, we also extracted and interpreted simple yet important network motifs that involve
    the interplays between TF-miRNA and their targets. In summary, RegNetwork provides an
    integrated resource on the prior information for gene regulatory relationships, and it enables
    us to further investigate context-specific transcriptional and post-transcriptional regulatory
    interactions based on domain-specific experimental data.

    '''
    def __init__(self,
                 species='hsa',
                 directions=True,
                 sources=False,
                 exclude=None,
                 include=None,
                 make_edges_unique=False,
                 root_url=None,
                 local_path=None,
                 verbose=True,
                 node_columns=(1,3),
                 annotate=True,
                 databases_as_list=True):
        '''
        Args:

            species (str): Species for which you want the RegNetwork graph for. Available
                           are 'hsa' (human) and 'mmu' (mouse).
                           Default: 'hsa'
            directions (bool): Whether to include the direction of the edges where known.
                               Default: True
            sources (bool): ...
            exclude (list): List of edge types to exclude.
                            Default: []
            include (list): List of edge types to include.
                            Default: []
            make_edges_unique (bool): If True edges with identical incident nodes will result
                                      in only one edge and the attributes will be accessible
                                      as a list. If False, the graph can contain multiedges. 
                                      For the DeRegNet algorithms any is fine, it will mostly
                                      be a matter of downstream convenience whether you choose
                                      one or the other.
                                      Default: False
            root_url (str): Base URL from where you can access the RegNetwork data for download
                            Default: 'http://www.regnetworkweb.org/download'
            local_path (str): Path where to store downloaded files.
                              Default: '${HOME}/.deregnet/graphs/regnetwork'
            verbose (bool): If True you get some additional messages during download, etc.
            node_columns (tuple): DO NOT USE THIS ARGUMENT, THERE REALLY SHOULD BE NO REASON
            annotate (bool): If True the graph will have additional attributes, like several
                             ID systems identifiers as node attribute, etc.
                             Default: True
            databases_as_list (bool): If True the database origin attribute of the edges will be
                                      in list format. Otherwise it will be a comma-seperated string.
                                      Default: True
        '''
        root_url = DEFAULT_REG_NETWORK_DOWNLOAD_URL if root_url is None else root_url
        local_path = DEFAULT_REG_NETWORK_LOCAL_PATH if local_path is None else local_path
        data = self.get(species, directions, sources, local_path, root_url=root_url, verbose=verbose)
        attributes = {}
        if sources and directions:
            attributes[4] = 'direction' if not make_edges_unique else 'directions'
        elif not sources:
            attributes[4] = 'databases'
            attributes[5] = 'evidence' if not make_edges_unique else 'evidences'
            attributes[6] = 'confidence' if not make_edges_unique else 'confidences'
        source_column, target_column = node_columns
        igraph_init_kwargs = table_to_igraph_init_kwargs(data,
                                                         source_column=source_column,
                                                         target_column=target_column,
                                                         attributes=attributes,
                                                         directed=True,
                                                         exclude=exclude,
                                                         include=include,
                                                         make_edges_unique=make_edges_unique)
        super().__init__(make_edges_unique, **igraph_init_kwargs)
        self.species = species
        self.directions = directions
        if annotate and node_columns == (1,3):
            self.annotate(databases_as_list)

    def annotate(self, databases_as_list):
        self.vs['mirna'] = [True if v['name'].startswith('MI') else False for v in self.vs]
        self.vs['protein_coding'] = [not v for v in self.vs['mirna']]
        tfs = {e.source for e in self.es}
        self.vs['tf'] = [(v.index in tfs) if v['protein_coding'] else False for v in self.vs]
        self.vs['node_type'] = ['gene' if v['protein_coding'] else 'mirna' for v in self.vs]
        self.vs['node_type'] = [self.vs['node_type'][i] if not v['tf'] else 'tf' for i, v in enumerate(self.vs)]
        self.map_nodes(BioMap().get_mapper('mirbase'), TO='alias', target_attr='mirna_alias')
        # TODO: map MiRBase families
        if self.species == 'hsa':
            self.map_nodes_to_multiple_targets(BioMap().get_mapper('hgnc'),
                                               FROM='entrez_id',
                                               TO=['entrez_id', 'ensembl_id', 'symbol'],
                                               target_attrs=['entrez', 'ensembl', 'symbol'])
        else:
            self.map_nodes_to_multiple_targets(BioMap().get_mapper('mgi_entrez'),
                                               FROM='entrez',
                                               TO=['entrez', 'symbol', 'name'],
                                               target_attrs=['entrez', 'symbol', 'name'])
            self.map_nodes(BioMap().get_mapper('mgi_ensembl'),
                           FROM='symbol', TO='ensembl',
                           source_attr='symbol', target_attr='ensembl')
            self.vs['name'] = [','.join(v['mirna_alias']).split(',')[0].split(self.species+'-')[-1] if v['mirna'] else v['symbol'] for v in self.vs]
        # edges
        if databases_as_list:
            self.es['databases'] = [s.split(',') for s in self.es['databases']]
        if self.directions:
            g = RegNetwork(species=self.species, directions=True, sources=True, annotate=False)
            edges = {(g.vs[e.source]['name'], g.vs[e.target]['name']): e['direction'] for e in g.es}
            self.es['direction'] = ['-/-' if (self.vs[e.source]['name'], self.vs[e.target]['name']) not in edges
                                    else edges[(self.vs[e.source]['name'], self.vs[e.target]['name'])]
                                    for e in self.es]
            self.es['direction'] = ['--|' if self.vs[e.source]['mirna'] else e['direction'] for e in self.es]
            self.annotate_with_edge_types()

    def annotate_with_edge_types(self):
        def get_edge_type(self, edge):
            source = self.vs[edge.source]
            target = self.vs[edge.target]
            if source['mirna']:
                if target['tf']:
                    edge_type = 'mirna-tf'
                elif target['mirna']:
                    edge_type = 'mirna-mirna'
                else:
                    edge_type = 'mirna-gene'
            elif source['tf']:
                if target['tf']:
                    edge_type = 'tf-tf'
                elif target['mirna']:
                    edge_type = 'tf-mirna'
                else:
                    edge_type = 'tf-gene'
            return edge_type

        self.es['edge_type'] = [get_edge_type(self, edge) for edge in self.es]

    @classmethod
    def download(cls, what='RegulatoryDirections', verbose=True, root_url=None, local_path=None):
        root_url = DEFAULT_REG_NETWORK_DOWNLOAD_URL if root_url is None else root_url
        if what == 'RegulatoryDirections':
            filename = what+'.rar'
        else:
            filename = what+'.zip'
        url = root_url+'/'+filename
        local_path = DEFAULT_REG_NETWORK_LOCAL_PATH if local_path is None else local_path
        local_file = os.path.join(local_path, filename)
        cls._download(url, local_file, verbose)

    @classmethod
    def get(cls, species='hsa', directed=True, sources=False, local_path=None, **kwargs):
        # ---
        import rarfile
        # ---
        species = 'human' if species == 'hsa' else 'mouse'
        what = species if not directed else 'RegulatoryDirections'
        local_path = DEFAULT_REG_NETWORK_LOCAL_PATH if local_path is None else local_path
        if directed:
            local_file = os.path.join(local_path, what+'.rar')
        else:
            local_file = os.path.join(local_path, species+'.zip')
        if not sources:
            return cls.get_data_from_form(species, local_path, skiprows=1, header=None)
        if not os.path.isfile(local_file):
            cls.download(what, local_path=local_path, **kwargs)
        if directed:
            filename = 'kegg.'+species+'.reg.direction'
            with rarfile.RarFile(local_file) as rf:
                with rf.open(filename) as fp:
                    return pd.read_table(fp, sep='\s+', skiprows=1, header=None)
        else:
            filename = species+'.source'
            with zipfile.ZipFile(local_file) as zf:
                with zf.open(filename) as fp:
                    return pd.read_table(fp, header=None, low_memory=False)

    @classmethod
    def get_data_from_form(cls, species='hsa', local_path=None, **kwargs):
        local_path = DEFAULT_REG_NETWORK_LOCAL_PATH if local_path is None else local_path
        local_file = cls.download_via_form(species, local_path)
        return pd.read_csv(local_file, low_memory=False, **kwargs)

    @classmethod
    def download_via_form(cls, species='hsa', local_path=None):
        def response_status(response, response_name):
            if response.status_code != 200:
                print('ERROR during %s request!' % response_name)
                return None

        local_path = DEFAULT_REG_NETWORK_LOCAL_PATH if local_path is None else local_path
        species = 'human' if species in {'hsa', 'human'} else 'mouse'
        local_file = os.path.join(local_path, species+'.csv')
        if os.path.isfile(local_file):
            return local_file
        if species == 'mouse':
            search_response = requests.get(cls.search_request('mmu'))
            if response_status(search_response, 'search'): return None
            export_response = requests.get(cls.export_request('mmu'))
            if response_status(search_response, 'export'): return None
            with open(local_file, 'wb') as fp:
                fp.write(export_response.content)
        else:
            search_response = requests.get(cls.search_request('hsa', 'Experimental'))
            if response_status(search_response, 'search'): return None
            export_response = requests.get(cls.export_request('hsa'))
            if response_status(search_response, 'export'): return None
            local_file_experimental = os.path.join(local_path, 'human_experimental.csv')
            with open(local_file_experimental, 'wb') as fp:
                fp.write(export_response.content)
            search_response = requests.get(cls.search_request('hsa', 'Predicted'))
            if response_status(search_response, 'search'): return None
            export_response = requests.get(cls.export_request('hsa'))
            if response_status(search_response, 'export'): return None
            local_file_predicted = os.path.join(local_path, 'human_predicted.csv')
            with open(local_file_predicted, 'wb') as fp:
                fp.write(export_response.content)
            experimental = pd.read_csv(local_file_experimental, low_memory=False)
            predicted = pd.read_csv(local_file_predicted, low_memory=False)
            human = pd.concat([experimental, predicted])
            human.to_csv(local_file, index=False)
        return local_file

    @classmethod
    def search_request(cls, species='hsa', evidence='all'):
        url = 'http://www.regnetworkweb.org/search.jsp?'
        url += 'searchItem=&searchType=all&'
        url += 'organism='+('human&' if species == 'hsa' else 'mouse&')
        url += 'database=all&evidence='+evidence+'&confidence=all&'
        url += 'resultsPerPage=30&prevValidPN=1&orderBy=RegSymbol_Asc&pageNumber=1'
        return url

    @classmethod
    def export_request(cls, species='hsa'):
        url = 'http://www.regnetworkweb.org/export.jsp?format=csv&'
        url += 'sql=SELECT+*+FROM+'+('human' if species == 'hsa' else 'mouse')+'+WHERE'
        url += '+%271%27+ORDER+BY+UPPER%28regulator_symbol%29+ASC'
        return url
