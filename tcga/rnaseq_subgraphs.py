'''
Finds deregulated subgraphs based on RNASeq gene expression in TCGA
'''

import os
import argparse

import igraph as ig

from biomap import BioMap
from deregnet.core import SubgraphFinder
from deregnet_tcga.scores import (get_rnaseq_score,
                                  get_global_rnaseq_score)
from deregnet_tcga.layers import Layers
from deregnet_tcga.graph import (get_expression_induced_subgraph,
                                 prepare_expression_indicator)


__FILEDIR__ = os.path.dirname(os.path.abspath(__file__))

GRAPH_PATH = os.path.join(__FILEDIR__, 'graph', 'kegg_hsa.graphml')

PATIENT_SPECIFIC_LAYERS = ['snv_as_roots', 'snv_as_terminals']


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-l', '--layer', type=str, default='null',
                        help='Which predefined layer to use')
    parser.add_argument('-d', '--dataset', type=str, default='uvm',
                        help='Which TCGA dataset to use')
    parser.add_argument('-m', '--mode', type=str, default='deregulated',
                        help='Which mode to choose: deregulated, upregulated,'
                        ' downregulated')
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
    parser.add_argument('-n', '--normalize-wrt', type=str, default='control',
                        help='Normalization strategy for RNASeq data')
    parser.add_argument('--suboptimal', action='store_true',
                        help='Whether to find next best suboptimal subgraph.')
    parser.add_argument('--max-suboptimal', type=int, default=0,
                        help='Whether to find next best suboptimal subgraph.')
    parser.add_argument('--global', dest='_global', action='store_true',
                        help='Find global subgraphs.')
    parser.add_argument('--global-score', dest='global_score', type=str, default='log2fold',
                        help='Log2Fold or Binary.')
    return parser.parse_args()


def prepare_rnaseq_score(dataset, compare_to):
    id_mapper = BioMap().get_mapper('hgnc')
    rnaseq_score = get_rnaseq_score(dataset, compare_to=compare_to)
    rnaseq_score.index = [gene.split('.')[0] for gene in rnaseq_score.index]
    rnaseq_score.index = list(id_mapper.map(list(rnaseq_score.index),
                                            FROM='ensembl',
                                            TO='entrez'))
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


def read_fails(base_path):
    try:
        open(os.path.join(base_path, 'failed.txt'), 'x')
    except:
        pass
    with open(os.path.join(base_path, 'failed.txt'), 'r') as fp:
        fails = [fail.strip() for fail in fp.read().split('\n')]
        return [fail for fail in fails if fail]
    return []

def global_rnaseq_subgraphs(args):
    graph = ig.Graph.Read_GraphML(GRAPH_PATH)
    mode, orientation = args.mode.split('/')
    flip = True if orientation == 'upstream' else False
    base_path = os.path.join('rnaseq/global',
                             args.dataset,
                             args.global_score,
                             mode,
                             orientation,
                             args.layer)
    receptors, terminals = Layers.get(args.layer)
    if not os.path.isdir(base_path):
        os.makedirs(base_path)
    score = get_global_rnaseq_score(args.dataset, args.global_score)
    abs_vals, minmax = get_mode_args(mode)
    finder = SubgraphFinder(graph)
    try:
        result = finder.run_average_deregnet(score,
                                             min_size=args.min_size,
                                             max_size=50,
                                             time_limit=args.time_limit,
                                             gap_cut=args.gap_cut,
                                             abs_values=abs_vals,
                                             model_sense=minmax,
                                             flip_orientation=flip,
                                             receptors=receptors,
                                             terminals=terminals,
                                             num_suboptimal=args.max_suboptimal,
                                             max_overlap=0.2,
                                             min_num_terminals=args.min_num_terminals)
    except:
        print('could not find subgraphs')
        return
    try:
        result.to_graphml(path=base_path)
    except:
        print('could not write subgraphs.')
        return

def local_rnaseq_subgraphs(args):
    hgnc = BioMap().get_mapper('hgnc')
    graph = ig.Graph.Read_GraphML(GRAPH_PATH)
    #
    rnaseq_score = prepare_rnaseq_score(args.dataset,
                                        compare_to=args.normalize_wrt)
    expression_ind_tumor = prepare_expression_indicator(args.dataset, threshold=100, tumor=True)
    expression_ind_control = prepare_expression_indicator(args.dataset, threshold=100, tumor=False)
    patients = list(rnaseq_score.columns)
    abs_vals, minmax = get_mode_args(args.mode)
    base_path = os.path.join('rnaseq', args.layer, args.mode, args.dataset)
    if not os.path.isdir(base_path):
        os.makedirs(base_path)
    failed = set(read_fails(base_path))
    for patient_id in patients:
        path = os.path.join(base_path, patient_id)
        if not os.path.isdir(path):
            os.makedirs(path)
        if patient_id in failed:
            print('Ignoring previously failed patient %s.' % patient_id)
            continue
        suboptimality_index = 0
        suboptimality_excludes = set()

        print(path)
        print(os.listdir(path))
        subgraphs = os.listdir(path)
        suboptimal_subgraphs = set(subgraphs) - {'optimal.graphml'}
        if not os.path.isfile(os.path.join(path, 'optimal.graphml')):
            suboptimality_index = 0
        else:
            print({f.split('_')[1].split('.')[0] for f in suboptimal_subgraphs})
            suboptimality_index = max({int(f.split('_')[1].split('.')[0]) for f in suboptimal_subgraphs} | {0}) + 1
        print('Suboptimality index: %d' % suboptimality_index)
        for f in subgraphs:
            fpath = os.path.join(path, f)
            nodes = set(ig.Graph.Read_GraphML(fpath).vs['entrez'])
            suboptimality_excludes = suboptimality_excludes.union(nodes)
        if (suboptimality_index == 1 and not args.suboptimal):
            print('optimal found. Not going for suboptimal_1.')
            continue
        if (suboptimality_index > 1 and not args.suboptimal):
            print('suboptimal_%d found. Not going for suboptimal_%d.' %(suboptimality_index-1, suboptimality_index))
            continue
        if (suboptimality_index > args.max_suboptimal):
            print('Suboptimality index exceeding maximum suboptimality index.')
            continue
        try:
            expression_induced_graph = \
                get_expression_induced_subgraph(graph,
                                                patient_id,
                                                expression_ind_tumor,
                                                expression_ind_control)
            print(patient_id,
                  ': #nodes ',
                  str(len(expression_induced_graph.vs)),
                  ' #edges ',
                  str(len(expression_induced_graph.es)))
        except Exception:
            failed.add(patient_id)
            print('No mRNA expression for ', patient_id)
            continue
        finder = SubgraphFinder(expression_induced_graph)
        score = rnaseq_score[patient_id].to_dict()
        if args.layer in PATIENT_SPECIFIC_LAYERS:
            receptors, terminals = Layers.get(args.layer,
                                              cancer_type=args.dataset.upper(),
                                              patient=patient_id)
        else:
            receptors, terminals = Layers.get(args.layer)
        flip = True if args.layer.startswith('terminal') else False
        if args.layer == 'snv_as_terminals':
            flip = True
        if args.layer.startswith('terminal') or args.layer.startswith('rooted'):
            if receptors and len([g for g in receptors if expression_induced_graph.vs.select(entrez_eq=g)]) == 0:
                continue
            if terminals and len([g for g in terminals if expression_induced_graph.vs.select(entrez_eq=g)]) == 0:
                continue
        result = None
        exclude = hgnc.map(args.exclude.split(','), TO='entrez')
        exclude = set(exclude).union(suboptimality_excludes)
        if args.layer in {'snv_as_terminals', 'snv_as_roots'}:
            receptor_excludes = set(receptors) if receptors else set()
            terminal_excludes = set(terminals) if terminals else set()
            exclude = exclude.difference(set(receptor_excludes), set(terminal_excludes))
        exclude = [gene for gene in exclude if gene]
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
                                                 excluded_nodes=exclude)
        except:
            failed.add(patient_id)
        try:
            name = 'optimal.graphml'
            if suboptimality_index > 0:
                name = 'suboptimal_'+str(suboptimality_index)+'.graphml'
            result.optimal.write_graphml(os.path.join(path, name))
        except:
            failed.add(patient_id)
    write_fails(failed, base_path)


def main(args):
    if args._global:
        global_rnaseq_subgraphs(args)
    else:
        local_rnaseq_subgraphs(args)

if __name__ == '__main__':
    main(parse_args())
