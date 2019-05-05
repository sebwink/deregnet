import os
import numpy as np
import pandas as pd
import igraph as ig
import matplotlib.pyplot as plt
import seaborn as sns

DATA_HOME = os.environ.get('DEREGNET_TCGA_DATA_HOME',
                           os.path.join(os.environ['HOME'], 'projects/DeRegNet/deregnet/tcga'))
BASE_GRAPH = os.path.join(DATA_HOME, 'graph/kegg_hsa.graphml')

class PatientSubgraphs:
    def __init__(self, path=None, base_graph=BASE_GRAPH):
        if isinstance(path, str):
            path = [path]
        if path is None:
            self._graphs = []
        else:
            self._graphs = self.read_graphs(path[0])
            for p in path[1:]:
                self._graphs = self._graphs + self.read_graphs(path)
        self._base_graph = base_graph

    @classmethod
    def from_graph_list(self, graphs, base_graph=BASE_GRAPH):
        patient_subgraphs = PatientSubgraphs(None, base_graph)
        patient_subgraphs._graphs = graphs
        return patient_subgraphs

    @classmethod
    def read_graphs(cls, path):
        graphs = { }
        for f in os.listdir(path):
            fpath = os.path.join(path, f)
            if not fpath.endswith('.graphml'):
                continue
            if f.startswith('optimal'):
                graphs[0] = ig.Graph.Read_GraphML(fpath)
            else:
                i = int(f.split('_')[1].split('.')[0])
                graphs[i] = ig.Graph.Read_GraphML(fpath)
        return [graphs[i] for i, _ in enumerate(graphs)]

    def __getitem__(self, i):
        if i >= len(self._graphs):
            raise IndexError
        return self._graphs[i]

    def __add__(self, other):
        assert self._base_graph == other._base_graph
        graphs = self._graphs + other._graphs
        return self.from_graph_list(graphs, self._base_graph)

    def nodes(self, i, attr='symbol', if_attr_true=None):
        if if_attr_true is None:
            return {v[attr] for v in self[i].vs}
        else:
            return {v[attr] for v in self[i].vs if v.attributes().get(if_attr_true, False)}

    def subgraph_score(self, i, mode='average', absval=False):
        if len(self._graphs) == 0:
            return 0.0
        if absval:
            score_sum = sum(abs(v['deregnet_score']) for v in self[i].vs)
        else:
            score_sum = sum(v['deregnet_score'] for v in self[i].vs)
        if mode == 'average':
            return score_sum / len(self[i].vs)
        return score_sum

    def aggregate_subgraph_scores(self, aggregation_mode='average', mode='average', absval=False):
        aggregated_score_sum = sum(self.subgraph_score(i, mode, absval) for i in range(len(self._graphs)))
        if aggregation_mode == 'average':
            return aggregated_score_sum / len(self._graphs)
        return aggregated_score_sum

    def node_union(self, attr='symbol', if_attr_true=None):
        if len(self._graphs) == 0:
            return set()
        return self.nodes(0, attr, if_attr_true).union(*[self.nodes(i, attr, if_attr_true) for i in range(1, len(self._graphs))])

    def node2attr(self, attr, node_id='symbol'):
        node2attr = []
        for subgraph in self._graphs:
            node2attr.append({v[node_id]: v[attr] for v in subgraph.vs})
        ret = { }
        for n2a in node2attr:
            ret = { **ret, **n2a }
        return ret

    def union_graph(self,
                    base_graph=None,
                    match_attr='symbol'):
        if base_graph is None:
            base_graph = self._base_graph
        if isinstance(base_graph, str):
            bsgrph = ig.Graph.Read_GraphML(base_graph)
        else:
            bsgrph = base_graph
        nodes = eval('bsgrph.vs.select('+match_attr+'_in=self.node_union(match_attr))')
        subgraph = bsgrph.subgraph(nodes)
        if len(self._graphs) > 0:
            for attr in self._graphs[0].vs.attribute_names():
                subgraph.vs[attr] = [self.node2attr(attr, match_attr)[v[match_attr]] for v in subgraph.vs]
        return subgraph


class CohortSubgraphs:

    def __init__(self,
                 patient2subgraph,
                 base_graph=BASE_GRAPH):
        self._base_graph = base_graph
        self._patient2subgraph = patient2subgraph
        self._path = None

    def from_patient_subset(self, patients):
        patient2subgraph = {p: sg for p, sg in self._patient2subgraph.items() if p in patients}
        return CohortSubgraphs(patient2subgraph, self._base_graph)
    
    @classmethod
    def from_path_with_score_filter(cls, path, subset=None, threshold=0.3, base_graph=BASE_GRAPH, data_home=DATA_HOME):
                      
        def cmp_downreg(score, threshold=-0.3):
            return score < threshold

        def cmp_upreg(score, threshold=0.3):
            return score > threshold

        path = os.path.join(DATA_HOME, path)
        subgraphs = CohortSubgraphs.from_cohort_path(path, subset=subset, base_graph=base_graph, data_home=data_home)
        if 'downregulated' in path:
            cmp = cmp_downreg
        else:
            cmp = cmp_upreg
        deregulated = [patient for patient, score in subgraphs.subgraph_scores().items() if cmp(score)]
        subgraphs = CohortSubgraphs.from_cohort_path(path, subset=deregulated, base_graph=base_graph, data_home=data_home)
        return subgraphs
        
    @classmethod
    def from_cohort_path(cls,
                         cohort_path,
                         subset=None,
                         base_graph=BASE_GRAPH,
                         data_home=DATA_HOME):

        
        cohort_data = os.path.join(DATA_HOME, cohort_path)
        patient2subgraph = { }
        for patient_id in subset or os.listdir(cohort_data):
            path = os.path.join(cohort_data, patient_id)
            if os.path.isfile(path):
                continue
            patient2subgraph[patient_id] = PatientSubgraphs(path)
        subgraphs = CohortSubgraphs(patient2subgraph, base_graph)
        subgraphs._path = cohort_data
        return subgraphs

    def __add__(self, other):
        assert self._base_graph == other._base_graph
        patients = set(self.patients) | set(other.patients)
        patient2subgraphs = {
            p: self._patient2subgraph.get(p, PatientSubgraphs.from_graph_list([], self._base_graph)) + \
               other._patient2subgraph.get(p, PatientSubgraphs.from_graph_list([], self._base_graph))
            for p in patients
        }
        return CohortSubgraphs(patient2subgraphs, self._base_graph)

    def __or__(self, other):
        assert self._base_graph == other._base_graph
        return CohortSubgraphs({
                                 ** self._patient2subgraph,
                                 ** other._patient2subgraph   # TODO: handle case where patients appear in both: self and other ...
                               },
                               self._base_graph)

    def __getitem__(self, patient_id):
        return self._patient2subgraph[patient_id]

    @property
    def patients(self):
        return list(self._patient2subgraph.keys())

    def patient2nodes(self, i=None, patients=None, attr='symbol', if_attr_true=None):
        patients = patients if patients else self.patients
        if i is None:
            return {p: self[p].node_union(attr, if_attr_true) for p in patients }
        else:
            return {p: self[p].nodes(i, attr, if_attr_true) for p in patients}

    def node_histogram(self, i=None, patients=None, attr='symbol'):
        patients2nodes = self.patient2nodes(i, patients, attr)
        node_union = set().union(*[set(nodes) for nodes in patients2nodes.values()])
        histogram = {node: 0 for node in node_union}
        for p, nodes in patients2nodes.items():
            for node in nodes:
                histogram[node] += 1
        histogram = pd.DataFrame({'gene': list(histogram.keys()), 'count': list(histogram.values())})
        histogram.sort_values(by=['count'], ascending=False, inplace=True)
        return histogram

    def plot_node_histogram(self, patients=None, cutoff=5, attr='symbol', i=None):
        histogram = self.node_histogram(i, patients, attr)
        _, ax = plt.subplots(figsize=(14,8))
        g = sns.barplot(x='gene', y='count', data=histogram[histogram['count']>=cutoff], ax=ax)
        g.set_xticklabels(g.get_xticklabels(), rotation=90)

    def node_histogram_restricted_to_boolean_attr(self, rattr, i=None, patients=None, attr='symbol'):
        patients2nodes = self.patient2nodes(i, patients, attr, rattr)
        node_union = set().union(*[set(nodes) for nodes in patients2nodes.values()])
        histogram = {node: 0 for node in node_union}
        for p, nodes in patients2nodes.items():
            for node in nodes:
                histogram[node] += 1
        histogram = pd.DataFrame({'gene': list(histogram.keys()), 'count': list(histogram.values())})
        histogram.sort_values(by=['count'], ascending=False, inplace=True)
        return histogram

    def plot_node_histogram_restricted_to_boolean_attr(self, rattr, cutoff=5, **kwargs):
        histogram = self.node_histogram_restricted_to_boolean_attr(rattr, **kwargs)
        _, ax = plt.subplots(figsize=(14,8))
        g = sns.barplot(x='gene', y='count', data=histogram[histogram['count']>=cutoff], ax=ax)
        g.set_xticklabels(g.get_xticklabels(), rotation=90)


    def get_gene_cooccurence_matrix(self, i=None, patients=None, attr='symbol'):
        patient2nodes = self.patient2nodes(i, patients, attr)
        nodes = list({v for nodes in patient2nodes.values() for v in nodes})
        npairs = [(u,v) for i, u in enumerate(nodes) for v in nodes[i:]]
        cooccurence = {npair: 0 for npair in npairs}
        for patient, pnodes in patient2nodes.items():
            for n1, n2 in cooccurence:
                if n1 in pnodes and n2 in pnodes:
                    cooccurence[(n1,n2)] = cooccurence[(n1,n2)] + 1
        matrix = np.zeros((len(nodes), len(nodes)))
        for i, n1 in enumerate(nodes):
            for j, n2 in enumerate(nodes[i:]):
                val = cooccurence[(n1,n2)]
                matrix[i][j], matrix[j][i] = val, val
        return pd.DataFrame(data=matrix, index=nodes, columns=nodes)


    def subgraph_scores(self, i=0, mode='average', absval=None):
        if absval is None and self._path.split('/')[-2] == 'deregulated':
            absval = True
        else:
            absval = False
        return {patient : subgraphs.subgraph_score(i, mode, absval)
                for patient, subgraphs in self._patient2subgraph.items()}

    def union_graph(self, base_graph=None, attr='symbol'):
        if base_graph is None:
            base_graph = self._base_graph
        if isinstance(base_graph, str):
            bsgrph = ig.Graph.Read_GraphML(base_graph)
        else:
            bsgrph = base_graph
        return {patient_id: graphs.union_graph(bsgrph, attr)
                for patient_id, graphs in self._patient2subgraph.items()}

    def node_union(self, attr='symbol'):
        return {patient_id: graphs.node_union(attr)
                for patient_id, graphs in self._patient2subgraph.items()}

    def _matrix_setup(self, base_graph=None, subset=None, attr='symbol'):
        subset = self.patients if subset is None else subset
        p2g = self.union_graph(base_graph, attr)
        graphs, patients = list(p2g.values()), list(p2g.keys())
        n = len(graphs)
        return graphs, np.zeros((n,n)), patients

    def node_intersection_count_matrix(self, base_graph=None, scale=False, subset=None, attr='symbol'):
        graphs, matrix, patients = self._matrix_setup(base_graph, subset, attr)
        matrix = populate_matrix_symmetric(graphs, matrix, node_intersection_count, scale=scale)
        return pd.DataFrame(data=matrix, index=patients, columns=patients)

    def regulation_aware_intersection_count_matrix(self, base_graph=None, scale=False, subset=None, attr='symbol'):
        graphs, matrix, patients = self._matrix_setup(base_graph, subset, attr)
        matrix = populate_matrix_symmetric(graphs, matrix, regulation_aware_node_intersection_count, scale=scale)
        return pd.DataFrame(data=matrix, index=patients, columns=patients)

    def boolean_attr_overlap_count_matrix(self, attr, base_graph=None, scale=True, subset=None, match_attr='symbol'):
        graphs, matrix, patients = self._matrix_setup(base_graph, subset, match_attr)
        matrix = populate_matrix_symmetric(graphs, matrix, _boolean_attr_overlap, attr=attr, match_attr=match_attr, scale=scale)
        return pd.DataFrame(data=matrix, index=patients, columns=patients)

    def receptor_intersection_count_matrix(self, base_graph=None, scale=True, subset=None, attr='symbol'):
        return self.boolean_attr_overlap_count_matrix('deregnet_receptor', base_graph, scale, subset, attr)

    def terminal_intersection_count_matrix(self, base_graph=None, scale=True, subset=None, attr='symbol'):
        return self.boolean_attr_overlap_count_matrix('deregnet_terminal', base_graph, scale, subset, attr)
    
    def with_any_of(self, genes):
        genes = [genes] if isinstance(genes, str) else genes
        patient2nodes = self.patient2nodes()
        return list({patient for patient, nodes in patient2nodes.items()
                     if any([(gene in nodes) for gene in genes])})

    def with_none_of(self, genes):
        genes = [genes] if isinstance(genes, str) else genes
        patient2nodes = self.patient2nodes()
        return list({patient for patient, nodes in patient2nodes.items()
                     if not any([(gene in nodes) for gene in genes])})

    def with_all_of(self, genes):
        genes = [genes] if isinstance(genes, str) else genes
        patient2nodes = self.patient2nodes()
        return list({patient for patient, nodes in patient2nodes.items()
                     if all([(gene in nodes) for gene in genes])})

    def plot_node_histogram_for_subgraphs_with(self, genes, cutoff=30, _all=False):
        if _all:
            patients = get_subgraphs_with_all(genes)
        else:
            patients = get_subgraphs_with_any(genes)
        print('# patients: ', len(patients))
        self.plot_node_histogram(cutoff=cutoff, patients=patients)
    

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
    if len(nodes1) == 0 and len(nodes2) == 0:
        return 0.0
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

def regulation_aware_node_intersection_count(G1, G2, scale):
    nodes1 = { v['name'] for v in G1.vs }
    nodes2 = { v['name'] for v in G2.vs }
    if len(nodes1) == 0 or len(nodes2) == 0:
        return 0.0
    intersection = nodes1.intersection(nodes2)
    intersection1 = G1.vs.select(name_in=list(intersection))
    intersection2 = G2.vs.select(name_in=list(intersection))
    intersection_size = len(intersection)
    intersection = sum([v1['deregnet_score']*v2['deregnet_score'] for v1,v2 in zip(intersection1, intersection2)])
    if scale:
        ret = intersection / (len(nodes1) + len(nodes2) - intersection_size)
    else:
        ret = intersection
    return max(0, ret)

def regulation_aware_intersection_count_matrix(root, subset=None, scale=True):
    graphs, matrix, patients = read_graphs_setup_matrix(root, subset)
    return populate_matrix_symmetric(graphs, matrix, regulation_aware_node_intersection_count, scale=scale), patients

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

def _boolean_attr_overlap(G1, G2, attr, scale, match_attr):
    pos1 = {node[match_attr] for node in G1.vs if node[attr]}
    pos2 = {node[match_attr] for node in G2.vs if node[attr]}
    nodes1 = {node[match_attr] for node in G1.vs}
    nodes2 = {node[match_attr] for node in G2.vs}
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

def get_vogelstein(path='/home/sewi/projects/DeRegNet/deregnet/data/cancer_genes_vogelstein.txt'):
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
