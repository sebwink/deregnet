import os
import time
import shutil
import subprocess

import igraph as ig
import pandas as pd


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
        args = []
        if species != 'hsa':
            args = [species]
        mapper = id_mapper(*args)
        ids = mapper.map(ids, score_id_type, graph_id_type)
    df['ids'] = ids
    df = df[['ids', score_col]]
    df = df[df['ids'] != '']
    # df = df[str(df[score_col]) != '']
    return df 

def read_gmt(path2gmt):
    gmt = open(path2gmt, "r")
    geneSets = {}
    descriptions = {}
    for line in gmt:
        getthis = line.split("\n")[0]
        geneSet = getthis.split("\t")[0]
        description = getthis.split("\t")[1]
        genes = getthis.split("\t")[2:-1]
        genes = { gene for gene in genes if len(gene) > 0 }
        geneSets[geneSet] = genes
        descriptions[geneSet] = description
    gmt.close()
    return geneSets, descriptions

def parse_layer(layer_file, layer_from_gmt, graph_id_type, layer_id_type, species, id_mapper):
    layer = set(pd.read_table(layer_file)[[0]].tolist())
    if layer_from_gmt is not None:
        gmt_data = layer_from_gmt.split(',')
        gmt_file, gene_sets = gmt_data[0], gmt_data[1:]
        gmt, _ = read_gmt(gmt_file)
        for gene_set in gene_sets:
            layer |= set(gmt[gene_set])
            layer = list(layer)
        if layer_id_type != graph_id_type:
            args = []
            if species != 'hsa':
                args = [species]
            mapper = id_mapper(*args)
            layer = mapper.map(layer, layer_id_type, graph_id_type)
    return layer

def write_tmp_files(graph,
                    graph_id_attr,
                    graph_id_type,
                    path2score,
                    score_id_type,
                    sep,
                    id_col,
                    score_col,
                    species,
                    id_mapper,
                    receptors_file,
                    receptors_from_gmt,
                    receptor_id_type,
                    terminals_file,
                    terminals_from_gmt,
                    terminal_id_type,
                    **kwargs):
    '''
    Write all necessary temporary files.
    '''
    tmp_files = TmpFileHandle()
    # graph 
    graph = remove_self_loops(graph)
    tmp_files.write_graph(graph, graph_id_attr)
    # score
    scores = parse_scores(path2score,
                          id_col,
                          score_col,
                          sep,
                          score_id_type,
                          graph_id_type,
                          species,
                          id_mapper,
                          **kwargs)
    tmp_files.write_scores(scores)
    # receptors
    receptors = parse_layer(receptors_file,
                            receptors_from_gmt,
                            graph_id_type,
                            receptor_id_type,
                            species,
                            id_mapper)
    tmp_files.write_receptors(receptors)
    # terminals
    terminals = parse_layer(terminals_file,
                            terminals_from_gmt,
                            graph_id_type,
                            terminal_id_type,
                            species,
                            id_mapper)
    tmp_files.write_terminals(terminals)
    # 
    return tmp_files

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
            return (':'.join(first_line.split(':')[1:])).strip()

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

#####################################################################################

#####################################################################################

def time_stamp():
    return time.strftime('%Y-%m-%d-%H-%M-%S', time.gmtime())

DEREGNET_TMPPATH = os.path.join(os.path.expanduser('~'), '.deregnet/tmp')

DEREGNET_BINPATH=os.path.join(__file__,'../../bin')

def igraph_to_lgf(graph,
                  path,
                  additional_graph_attributes = {},
                  id_attr = 'id'):
    '''
    Writes an igraph graph to lemon graph format (lgf)
    (discarding all node and edge attributes).
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


class AbsoluteDeregnetArguments:
    pass


class AverageDeregnetArguments:
    pass


class SubgraphFinder:
    '''

    '''
    def __init__(self, graph,
                       id_attr='name',
                       deregnet_binpath=None,
                       tmp_file_path=DEREGNET_TMPPATH):
        '''

        '''
        self.deregnet_binpath = DEREGNET_BINPATH if deregnet_binpath is None else deregnet_binpath
        tmp_file_path = DEREGNET_TMPPATH if tmp_file_path is None else tmp_file_path
        self.tmp_file_path = os.path.join(tmp_file_path, time_stamp())
        if not os.path.isdir(self.tmp_file_path):
            os.makedirs(self.tmp_file_path)
        self.graph_file = os.path.join(self.tmp_file_path, 'graph.lgf')
        self.graph = graph
        self.id_attr = id_attr
        self._graph_to_lgf()  # writes self.graph_file

    def run(self, args):
        '''

        '''
        if isinstance(args, AverageDeregnetArguments):
            return self.run_average_deregnet(**args())
        elif isinstance(arg, AbsoluteDeregnetArguments):
            return self.run_absolute_deregnet(**args())
        return None

    def run_average_deregnet(self,
                             scores={},
                             default_score=None,
                             receptors=None,
                             terminals=None,
                             excluded_nodes=None,
                             included_nodes=None,
                             flip_orientation=False,
                             min_size=15,
                             max_size=50,
                             num_suboptimal=0,
                             max_overlap=0,
                             abs_values=False,
                             model_sense='max',
                             algorithm='GeneralizedCharnesCooper',
                             time_limit=1200,
                             gap_cut=0.05,
                             debug=False):
        '''

        '''
        # set temporary paths and write temporary files
        rundir = os.path.join(self.tmp_file_path, 'run_'+time_stamp())
        os.makedirs(rundir)
        # temporary files
        score_file = os.path.join(rundir, 'scores.tsv')
        terminals_file = os.path.join(rundir, 'terminals.txt')
        receptors_file = os.path.join(rundir, 'receptors.txt')
        exclude_file = os.path.join(rundir, 'exclude.txt')
        include_file = os.path.join(rundir, 'include.txt')
        # handle scores
        self._prepare_scores(score_file, scores, default_score)  # writes self.score_file
        # handle algorithm synonyms
        if algorithm == 'GeneralizedCharnesCooper':
            algorithm = 'gcc'
        elif algorithm == 'Dinkelbach':
            algorithm = 'dta'
        elif algorithm == 'ObjectiveVariableTransform':
            algorithm = 'ovt'
        # assemble arguments
        avgdrgnt_argdict = {
                             '--graph' : self.graph_file,
                             '--score' : score_file,
                             '--output-dir' : rundir,
                             '--gap-cut' : str(gap_cut),
                             '--min-size' : str(min_size),
                             '--max-size': str(max_size),
                             '--suboptimal' : str(num_suboptimal),
                             '--max-overlap-percentage' : str(max_overlap),
                             '--time-limit' : str(time_limit),
                             '--model-sense' : model_sense,
                             '--algorithm' : algorithm
                            }
        # write receptor and terminal layers to temporary files
        if receptors:
            self._write_geneset(receptors_file, receptors)
            avgdrgnt_argdict['--receptors-file'] = self.receptors_file
        if terminals:
            self._write_geneset(terminals_file, terminals)
            avgdrgnt_argdict['--terminals-file'] = self.terminals_file
        # write temporary files for manually excluded or included nodes
        if excluded_nodes:
            self._write_geneset(exclude_file, excluded_nodes)
            avgdrgnt_argdict['--exclude-file'] = self.exclude_file
        if included_nodes:
            self._write_geneset(include_file, included_nodes)
            avgdrgnt_argdict['--receptors-file'] = self.receptors_file
        # arguments as list (for subprocess.call)
        avgdrgnt_args = []
        for arg in avgdrgnt_argdict:
            if avgdrgnt_argdict[arg] is not None:
                 avgdrgnt_args += [arg, avgdrgnt_argdict[arg]]
        # add flag arguments
        if flip_orientation:
            avgdrgnt_args += ['--flip-orientation']
        if abs_values:
            avgdrgnt_args += ['--absolute-values']
        # call C++ program
        avgdrgnt = os.path.join(self.deregnet_binpath, 'avgdrgnt')
        call = [avgdrgnt] + avgdrgnt_args
        if debug:
            print(call)
            subprocess.call(['gdb', '--args']+call)
        else:
            subprocess.call(call)  # TODO: reroute messages
        # parse back the results
        node_names_list = self._read_result(rundir)
        subgraphs = [self._get_subgraph(node_names) for node_names in node_names_list]
        self._write_deregnet_attrs(subgraphs, scores, default_score, receptors, terminals)
        return SubgraphFinderResult(optimal=subgraphs[0],
                                    suboptimal=subgraphs[1:],
                                    mode='avg')

    def run_absolute_deregnet(self, scores, root=None, flip_orientation=False, default_score=None):
        '''

        '''
        self._prepare_scores(scores, default_score)
        root = set() if root is None else {root}
        self._write_receptors(root)

    def _remove_self_loops(self, graph):
        self_loops = [e for e in graph.es if e.source == e.target]
        graph.delete_edges(self_loops)
        return graph

    def _graph_to_lgf(self):
        graph = self.graph.copy()
        graph = self._remove_self_loops(graph)
        igraph_to_lgf(graph, self.graph_file, id_attr=self.id_attr)

    def _prepare_scores(self, score_file, scores, default_score):
        if default_score is None:
            default_score = 0.0
        self._write_temporary_score_file(score_file, scores, default_score)

    def _write_temporary_score_file(self, score_file, scores, default_score):
        with open(score_file, 'w') as fp:
            for node in self.graph.vs:
                name = node[self.id_attr]
                fp.write(name+'\t'+str(float(scores.get(name, default_score)))+'\n')

    def _write_geneset(self, filepath, geneset):
        with open(filepath, 'w') as fp:
            for node in geneset:
                fp.write(node+'\n')

    def _read_result(self, rundir):
        output = os.path.join(rundir, 'plain')
        node_names_list = []
        for sif in os.listdir(output):
            node_names_list.append(self._nodenames_from_sif(os.path.join(output, sif)))
        return node_names_list

    def _nodenames_from_sif(self, sif):
        nodes = set()
        with open(sif, 'r') as sifp:
            for edge in sifp:
                edge = edge.split('\t')
                nodes.add(edge[0])
                nodes.add(edge[2].strip())  # there will be (!) no isolated nodes
        return nodes

    def _get_subgraph(self, node_names):
        nodes = self.graph.vs.select(name_in=node_names)
        return self.graph.subgraph(nodes, 'create_from_scratch')

    def _write_deregnet_attrs(self,
                              subgraphs,
                              score,
                              default_score,
                              receptors,
                              terminals):
        for subgraph in subgraphs:
            for node in subgraph:
                node['deregnet_score'] = self.scores.get(node['name'], default_score)
                node['deregnet_receptor'] = True if node['name'] in receptors else False
                node['deregnet_terminal'] = True if node['name'] in terminals else False

    def __del__(self):
        shutil.rmtree(self.tmp_file_path)


class SubgraphFinderResult:
    def __init__(self, optimal, suboptimal, mode):
        self.optimal = optimal
        self.suboptimal = suboptimal
        self._mode = mode

    @property
    def optimal_score(self):
        if self._mode == 'avg':
            return self.optimal_avg_score
        elif self._mode == 'abs':
            return self.optimal_abs_score

    @property
    def suboptimal_scores(self):
        if self.mode == 'avg':
            return self.suboptimal_avg_scores
        elif self.mode == 'abs':
            return self.suboptimal_abs_scores

    @property
    def scores(self):
        if self.mode == 'avg':
            return self.avg_scores
        elif self.mode == 'abs':
            return self.abs_scores

    @property
    def num_nodes_optimal(self):
        return len(self.optimal.vs)

    @property
    def num_nodes_suboptimal(self):
        return [len(subgraph.vs) for subgraph in self.suboptimal]

    @property
    def num_nodes(self):
        return [self.num_nodes_optimal] + self.num_nodes_suboptimal

    @classmethod
    def _num_nodes(cls, graph):
        return len(graph.vs)

    @classmethod
    def _abs_score(cls, graph):
        return sum(node['_deregnet_score'] for node in graph.vs)

    @classmethod
    def _avg_score(cls, graph):
        return cls._abs_score(graph) / cls._num_nodes(graph)

    @property
    def optimal_avg_score(self):
        return self._avg_score(self.optimal)

    @property
    def optimal_abs_score(self):
        return self._abs_score(self.optimal)

    @property
    def suboptimal_avg_scores(self):
        return [self._avg_score(subgraph) for subgraph in self.suboptimal]

    @property
    def suboptimal_abs_scores(self):
        return [self._abs_score(subgraph) for subgraph in self.suboptimal]

    @property
    def abs_scores(self):
        return [self.optimal_abs_score] + self.suboptimal_abs_scores

    @property
    def avg_scores(self):
        return [self.optimal_avg_score] + self.suboptimal_avg_scores

    @property
    def subgraphs(self):
        return [optimal] + suboptimal

    def to_graphml(path='.'):
        # TODO: stringify/ignore non-str,int,float attributes
        self.optimal.write_graphml(os.path.join(path, 'optimal.graphml'))
        for i, subgraph in enumerate(self.suboptimal):
            subgraph.write_graphml(os.path.join(path, 'suboptimal'+str(i+1)+'.graphml'))
