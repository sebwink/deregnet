'''

'''

import os
import time
import shutil
import subprocess

import igraph as ig
import pandas as pd

#####################################################################################

DEREGNET_TMPPATH = os.path.join(os.path.expanduser('~'), '.deregnet/tmp')

DEREGNET_BINPATH=os.path.join(os.path.dirname(__file__),'../../bin')

def time_stamp():
    return time.strftime('%Y-%m-%d-%H-%M-%S', time.gmtime())

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


def stringify_graph_attributes(_graph):

    def stringify(what):
        attrs = what.attributes()
        for item in what:
            for attr in attrs:
                val = item[attr]
                if not (isinstance(val, str) |
                        isinstance(val, int) |
                        isinstance(val, float)):
                    item[attr] = str(val)

    graph = _graph.copy()
    stringify(graph.vs)
    stringify(graph.es)
    return graph

class InvalidDeregnetArguments(Exception):
    pass

class DeregenetArguments:

    def __init__(self, scores={},
                       default_score=0.0,
                       excluded_nodes=None,
                       included_nodes=None,
                       flip_orientation=False,
                       num_suboptimal=0,
                       max_overlap=0.0,
                       abs_values=False,
                       model_sense='max',
                       time_limit=1200,
                       gap_cut=0.05,
                       debug=False):

        self.scores = dict(scores)
        self.default_score = float(default_score)
        self.excluded_nodes = excluded_nodes
        self.included_nodes = included_nodes
        self.flip_orientation = bool(flip_orientation)
        self.num_suboptimal = int(num_suboptimal)
        self._max_overlap = float(max_overlap)
        self.abs_values = bool(abs_values)
        self._model_sense = str(model_sense)
        self._time_limit = float(time_limit)
        self._gap_cut = float(gap_cut)
        self.debug = bool(debug)

        if not self.validate():
            raise InvalidDeregnetArguments

    def validate(self):
        valid = True
        if not self._validate_time_limit(self._time_limit):
            self.time_limit = self._time_limit
            valid = False
        if not self._validate_gap_cut(self._gap_cut):
            self.gap_cut = self._gap_cut
            valid = False
        if not self._validate_max_overlap(self._max_overlap):
            self.max_overlap = self._max_overlap
            valid = False
        if not self._validate_model_sense(self._model_sense):
            self.models_sense = self._model_sense
            valid = False
        return valid

    @property
    def max_overlap(self):
        return self._max_overlap

    @max_overlap.setter
    def max_overlap(self, percentage):
        if self._validate_max_overlap(percentage):
            self._max_overlap = percentage
        else:
            print('max_overlap must be in the interval [0,100]')

    def _validate_max_overlap(self, percentage):
        return 0 <= percentage <= 100

    @property
    def model_sense(self):
        return self._model_sense

    @model_sense.setter
    def model_sense(self, sense):
        if self._validate_model_sense(sense):
            self._model_sense = sense
        else:
            print('model_sense must be from {\'max\', \'min\'}')

    def _validate_model_sense(self, sense):
        return sense in {'min', 'max'}

    @property
    def time_limit(self):
        return self._time_limit

    @time_limit.setter
    def time_limit(self, limit):
        if self._validate_time_limit(limit):
            self._time_limit = limit
        else:
            print('time_limit must be None or a float >= 0')

    def _validate_time_limit(self, limit):
        return limit is None or limit >= 0

    @property
    def gap_cut(self):
        return self._gap_cut

    @gap_cut.setter
    def gap_cut(self, cut):
        if self._validate_gap_cut(cut):
            self._gap_cut = cut
        else:
            print('gap_cut must be None or a float from [0,1].')

    def _validate_gap_cut(self, cut):
        return cut is None or 0 < cut < 1

    def kwargs(self):
        return {
                 'scores': self.scores,
                 'default_score': self.default_score,
                 'excluded_nodes': self.excluded_nodes,
                 'included_nodes': self.included_nodes,
                 'flip_orientation': self.flip_orientation,
                 'num_suboptimal': self.num_suboptimal,
                 'max_overlap': self.max_overlap,
                 'abs_values': self.abs_values,
                 'model_sense': self.model_sense,
                 'time_limit': self.time_limit,
                 'gap_cut': self.gap_cut,
                 'debug': self.debug
               }


class AbsoluteDeregnetArguments(DeregenetArguments):

    def __init__(self, size=20, root=None, **kwargs):
        self._size = int(size)
        self._root = root
        super().__init__(**kwargs)

    def validate(self):
        valid = super().validate()
        if not self._validate_size(self._size):
            self.size = self._size
            valid = False
        if not self._validate_root(self._root):
            self.root = self._root
            valid = False
        return valid

    @property
    def size(self):
        return self._size

    @size.setter
    def size(self, size):
        if self._validate_size(size):
            self._size = size
        else:
            print('size should be >= 0')

    def _validate_size(self, size):
        return size >= 0

    @property
    def root(self):
        return self._root

    @root.setter
    def root(self, root):
        if self._validate_root(root):
            self._root = root
            if self.num_suboptimal > 0 and self.max_overlap < 1/self.size:
                self.max_overlap = 1/self.size
        else:
            print('root must be None or a string.')

    def _validate_root(self, root):
        return root is None or isinstance(root, str)

    def __call__(self):
        return {
                 **self.kwargs(),
                 **{
                     'size': self.size,
                     'root': self.root
                   }
               }

class AverageDeregnetArguments(DeregenetArguments):
    def __init__(self, min_size=15,
                       max_size=50,
                       receptors=None,
                       terminals=None,
                       algorithm='GeneralizedCharnesCooper',
                       **kwargs):
        self._min_size = int(min_size)
        self._max_size = int(max_size)
        self.receptors = receptors
        self.terminals = terminals
        self._algorithm = str(algorithm)
        super().__init__(**kwargs)

    def validate(self):
        valid = super().validate()
        if not self._validate_min_size(self._min_size):
            self.min_size = self._min_size
            valid = False
        if not self._validate_max_size(self._max_size):
            self.max_size = self._max_size
            valid = False
        if not self._validate_algorithm(self._algorithm):
            self.algorithm = self._algorithm
            valid = False
        return valid

    @property
    def max_size(self):
        return self._max_size

    @max_size.setter
    def max_size(self, size):
        if self._validate_max_size(size):
            self._max_size = size
        else:
            print('max_size must be None or an int >= min_size')

    def _validate_max_size(self, size):
        return size is None or size >= self.min_size

    @property
    def min_size(self):
        return self._min_size

    @min_size.setter
    def min_size(self, size):
        if self._validate_min_size(size):
            self._min_size = size
        else:
            print('min_size must be <= max_size')

    def _validate_min_size(self, size):
        return size <= self.max_size

    @property
    def algorithm(self):
        return self._algorithm

    @algorithm.setter
    def algorithm(self, name):
        if self._validate_algorithm(name):
            self._algorithm = name
        else:
            print('algorithm must be in {gcc, GeneralizedCharnesCooper, dta, Dinkelbach, ovt, ObjectiveVariableTransform}.')

    def _validate_algorithm(self, name):
        valid = name in {'gcc', 'GeneralizedCharnesCooper',
                         'dta', 'Dinkelbach',
                         'ovt', 'ObjectiveVariableTransform'}
        if self.gap_cut is not None and name in {'dta', 'Dinkelbach'}:
            print('It is not recommended to use the Dinkelbach-type algorithm together with a gap cut.'\
                  '\nHowever, for sake of curiosity, we allow it.')
        return valid

    def _validate_gap_cut(self, cut):
        valid = super()._validate_gap_cut(cut)
        if cut is not None and self.algorithm in {'dta', 'Dinkelbach'}:
            print('It is not recommended to use the Dinkelbach-type algorithm together with a gap cut.'\
                  '\nHowever, for sake of curiosity, we allow it.')
        return valid

    def __call__(self):
        return {
                 **self.kwargs(),
                 **{
                     'min_size': self.min_size,
                     'max_size': self.max_size,
                     'receptors': self.receptors,
                     'terminals': self.terminals,
                     'algorithm': self.algorithm
                   }
               }

class SubgraphFinder:
    '''

    '''
    def __init__(self, graph,
                       id_attr='name',
                       deregnet_binpath=None,
                       tmp_file_path=None,
                       delete_temporary_files=True):
        '''

        '''
        self.deregnet_binpath = DEREGNET_BINPATH if deregnet_binpath is None else deregnet_binpath
        tmp_file_path = DEREGNET_TMPPATH if tmp_file_path is None else tmp_file_path
        self.tmp_file_path = os.path.join(tmp_file_path, time_stamp())
        if not os.path.isdir(self.tmp_file_path):
            os.makedirs(self.tmp_file_path)
        self.delete_temporary_files = delete_temporary_files
        self.graph_file = os.path.join(self.tmp_file_path, 'graph.lgf')
        self.graph = graph
        self.id_attr = id_attr
        self._graph_to_lgf()  # writes self.graph_file

    def run(self, args):
        '''

        '''
        if isinstance(args, AverageDeregnetArguments):
            return self.run_average_deregnet(**(args()))
        elif isinstance(args, AbsoluteDeregnetArguments):
            return self.run_absolute_deregnet(**(args()))
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
                             time_limit=None,
                             gap_cut=None,
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
                             '--min-size' : str(min_size),
                             '--max-size': str(max_size),
                             '--suboptimal' : str(num_suboptimal),
                             '--max-overlap-percentage' : str(max_overlap),
                             '--model-sense' : model_sense,
                             '--algorithm' : algorithm
                            }
        if time_limit:
            avgdrgnt_argdict['--time-limit'] = str(time_limit)
        if gap_cut:
            avgdrgnt_argdict['--gap-cut'] = str(gap_cut)
        # write receptor and terminal layers to temporary files
        if receptors:
            self._write_geneset(receptors_file, receptors)
            avgdrgnt_argdict['--receptors-file'] = receptors_file
        if terminals:
            self._write_geneset(terminals_file, terminals)
            avgdrgnt_argdict['--terminals-file'] = terminals_file
        # write temporary files for manually excluded or included nodes
        if excluded_nodes:
            self._write_geneset(exclude_file, excluded_nodes)
            avgdrgnt_argdict['--exclude-file'] = exclude_file
        if included_nodes:
            self._write_geneset(include_file, included_nodes)
            avgdrgnt_argdict['--receptors-file'] = receptors_file
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
            print(' '.join(call))
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

    def run_absolute_deregnet(self,
                              scores={},
                              default_score=None,
                              root=None,
                              excluded_nodes=None,
                              included_nodes=None,
                              flip_orientation=False,
                              size=15,
                              num_suboptimal=0,
                              max_overlap=0,
                              abs_values=False,
                              model_sense='max',
                              time_limit=None,
                              gap_cut=None,
                              debug=False):
        # set temporary paths and write temporary files
        rundir = os.path.join(self.tmp_file_path, 'run_'+time_stamp())
        os.makedirs(rundir)
        # temporary files
        score_file = os.path.join(rundir, 'scores.tsv')
        exclude_file = os.path.join(rundir, 'exclude.txt')
        include_file = os.path.join(rundir, 'include.txt')
        # handle scores
        self._prepare_scores(score_file, scores, default_score)  # writes self.score_file
        # assemble arguments
        drgnt_argdict = {
                          '--graph' : self.graph_file,
                          '--score' : score_file,
                          '--output-dir' : rundir,
                          '--size' : str(size),
                          '--suboptimal' : str(num_suboptimal),
                          '--max-overlap-percentage' : str(max_overlap),
                          '--model-sense' : model_sense,
                        }
        if time_limit:
            drgnt_argdict['--time-limit'] = str(time_limit)
        if gap_cut:
            drgnt_argdict['--gap-cut'] = str(gap_cut)
        # write receptor and terminal layers to temporary files
        if root is not None:
            drgnt_argdict['--root'] = root
        # write temporary files for manually excluded or included nodes
        if excluded_nodes:
            self._write_geneset(exclude_file, excluded_nodes)
            drgnt_argdict['--exclude-file'] = exclude_file
        if included_nodes:
            self._write_geneset(include_file, included_nodes)
            drgnt_argdict['--include-file'] = receptors_file
        # arguments as list (for subprocess.call)
        drgnt_args = []
        for arg in drgnt_argdict:
            if drgnt_argdict[arg] is not None:
                 drgnt_args += [arg, drgnt_argdict[arg]]
        # add flag arguments
        if flip_orientation:
            drgnt_args += ['--flip-orientation']
        if abs_values:
            drgnt_args += ['--absolute-values']
        # call C++ program
        drgnt = os.path.join(self.deregnet_binpath, 'drgnt')
        call = [drgnt] + drgnt_args
        if debug:
            print(call)
            subprocess.call(['gdb', '--args']+call)
        else:
            subprocess.call(call)  # TODO: reroute messages
        # parse back the results
        node_names_list = self._read_result(rundir)
        subgraphs = [self._get_subgraph(node_names) for node_names in node_names_list]
        self._write_deregnet_attrs(subgraphs, scores, default_score, root=root) # TODO: register root
        return SubgraphFinderResult(optimal=subgraphs[0],
                                    suboptimal=subgraphs[1:],
                                    mode='abs')

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
        nodes = self.graph.vs.select(**{self.id_attr+'_in': node_names})
        return self.graph.subgraph(nodes, 'create_from_scratch')

    def _write_deregnet_attrs(self,
                              subgraphs,
                              scores,
                              default_score,
                              receptors=None,
                              terminals=None,
                              root=None):
        for subgraph in subgraphs:
            for node in subgraph.vs:
                node['deregnet_score'] = scores.get(node[self.id_attr], default_score)
                if receptors:
                    node['deregnet_receptor'] = True if node[self.id_attr] in receptors else False
                if terminals:
                    node['deregnet_terminal'] = True if node[self.id_attr] in terminals else False
                if root:
                    node['deregnet_root'] = True if node[self.id_attr] == root else False

    def __del__(self):
        if self.delete_temporary_files:
            shutil.rmtree(self.tmp_file_path)


class SubgraphFinderResult:
    def __init__(self, optimal, suboptimal, mode):
        self.optimal = optimal
        self.suboptimal = suboptimal
        self._mode = mode

    @property
    def mode(self):
        return self._mode

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
        return sum(node['deregnet_score'] for node in graph.vs)

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

    def to_graphml(self, path='.'):
        self.optimal.write_graphml(os.path.join(path, 'optimal.graphml'))
        for i, subgraph in enumerate(self.suboptimal):
            writable_subgraph = stringify_graph_attributes(subgraph)
            writable_subgraph.write_graphml(os.path.join(path, 'suboptimal'+str(i+1)+'.graphml'))
