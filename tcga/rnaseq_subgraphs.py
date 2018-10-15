import sys
import os
import argparse

import igraph as ig

from biomap import BioMap
from deregnet.core import SubgraphFinder
from deregnet_tcga.scores import get_rnaseq_score
from deregnet_tcga.layers import Layers
from deregnet_tcga.graph import get_expression_induced_subgraph, prepare_expression_indicator


__FILEDIR__ = os.path.dirname(os.path.abspath(__file__))

GRAPH_PATH = os.path.join(__FILEDIR__, 'graph', 'kegg_hsa.graphml')

PATIENT_SPECIFIC_LAYERS = ['genomic']

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-l', '--layer', type=str, default='null',
                        help='Which predefined layer to use')
    parser.add_argument('-d', '--dataset', type=str, default='uvm',
                        help='Which TCGA dataset to use')
    parser.add_argument('-m', '--mode', type=str, default='deregulated',
                        help='Which mode to choose: deregulated, upregulated, downregulated')
    parser.add_argument('-t', '--time-limit', type=int, default=600,
                        help='Maximal time to search for a subgraph')
    parser.add_argument('-g', '--gap-cut', type=float, default=None,
                        help='Gap cut to stop optimization prematurely')
    parser.add_argument('--min-size', type=int, default=10,
                        help='Minimal subgraph size.')
    parser.add_argument('--min-num-terminals', type=int, default=0,
                        help='Minimal number of terminals in subgraph.')
    parser.add_argument('--exclude', type=str, default='',
                        help='Which nodes to exclude from subgraphs.')
    return parser.parse_args()

def prepare_rnaseq_score(dataset):
    id_mapper = BioMap().get_mapper('hgnc')
    rnaseq_score = get_rnaseq_score(dataset, compare_to='control')
    rnaseq_score.index = [gene.split('.')[0] for gene in rnaseq_score.index]
    rnaseq_score.index = list(id_mapper.map(list(rnaseq_score.index), FROM='ensembl', TO='entrez'))
    return rnaseq_score

def get_mode_args(mode):
    if mode == 'deregulated':
        abs_vals = True
        minmax = 'max'
    else:
        abs_vals = False
        if mode == 'upregulated':
            minmax = 'max'
        elif mode == 'downregulated':
            minmax = 'min'
        else:
            raise ValueError
    return abs_vals, minmax

def write_fails(failed, base_path):
    with open(os.path.join(base_path, 'failed.txt'), 'w') as fp:
        for fail in failed:
            fp.write(fail+'\n')


def main(args):
    hgnc = BioMap().get_mapper('hgnc')
    graph = ig.Graph.Read_GraphML(GRAPH_PATH)
    #
    rnaseq_score = prepare_rnaseq_score(args.dataset)
    expression_ind = prepare_expression_indicator(args.dataset, threshold=100)
    patients = list(rnaseq_score.columns)
    abs_vals, minmax = get_mode_args(args.mode)
    base_path = os.path.join('rnaseq', args.layer, args.mode, args.dataset)
    if not os.path.isdir(base_path):
        os.makedirs(base_path)
    failed = []     # log patients for which no subgraph could be found
    for patient_id in patients:
        try:
            expression_induced_graph = get_expression_induced_subgraph(graph, patient_id, expression_ind)
            print(patient_id, ': #nodes ', str(len(expression_induced_graph.vs)), ' #edges ', str(len(expression_induced_graph.es)))
        except:
            failed.append(patient_id)
            print('No mRNA expression for ', patient_id)
        finder = SubgraphFinder(expression_induced_graph)
        path = os.path.join(base_path, patient_id)
        if not os.path.isdir(path):
            os.makedirs(path)
        if os.path.isfile(os.path.join(path, 'optimal.graphml')):
            continue
        score = rnaseq_score[patient_id].to_dict()
        if args.layer in PATIENT_SPECIFIC_LAYERS:
            receptors, terminals = Layers.get(args.layer, patient=patient_id)
        else:
            receptors, terminals = Layers.get(args.layer)
        flip = True if args.layer.startswith('terminal') else False
        if args.layer.startswith('terminal') or args.layer.startswith('rooted'):
            if receptors and len([g for g in receptors if expression_induced_graph.vs.select(entrez_eq=g)]) == 0:
                continue
            if terminals and len([g for g in terminals if expression_induced_graph.vs.select(entrez_eq=g)]) == 0:
                continue
        result = None
        print(args.exclude.split(','))
        exclude = hgnc.map(args.exclude.split(','), TO='entrez')
        print(exclude)
        try:
            result = finder.run_average_deregnet(score,
                                                 min_size=args.min_size,
                                                 max_size=50,
                                                 time_limit=args.time_limit,
                                                 gap_cut=args.gap_cut,
                                                 abs_values=abs_vals,
                                                 model_sense=minmax,
                                                 receptors=receptors,
                                                 terminals=terminals,
                                                 flip_orientation=flip,
                                                 min_num_terminals=args.min_num_terminals,
                                                 excluded_nodes=args.exclude.split(','))
        except:
            failed.append(patient_id)
        try:
            result.to_graphml(path)
        except:
            failed.append(patient_id)
    write_fails(failed, base_path)

if __name__ == '__main__':
    args = parse_args()
    main(args)
