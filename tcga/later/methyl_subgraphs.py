import sys
import os
import argparse

import igraph as ig

from biomap import BioMap
from deregnet.core import SubgraphFinder
from deregnet_tcga.scores import get_methylation_score
from deregnet_tcga.layers import Layers

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
    return parser.parse_args()

def prepare_methylation_score(dataset):
    id_mapper = BioMap().get_mapper('hgnc')
    rnaseq_score = get_methylation_score(dataset, normalize_wrt='genewise_median')
    #rnaseq_score.index = [gene.split('.')[0] for gene in rnaseq_score.index]
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
    graph = ig.Graph.Read_GraphML(GRAPH_PATH)
    finder = SubgraphFinder(graph)
    #
    methyl_score = prepare_methylation_score(args.dataset)
    patients = list(methyl_score.columns)
    abs_vals, minmax = get_mode_args(args.mode)
    base_path = os.path.join('methyl', args.layer, args.mode, args.dataset)
    if not os.path.isdir(base_path):
        os.makedirs(base_path)
    failed = []     # log patients for which no subgraph could be found
    for patient_id in patients:
        path = os.path.join(base_path, patient_id)
        if not os.path.isdir(path):
            os.makedirs(path)
        if os.path.isfile(os.path.join(path, 'optimal.graphml')):
            continue
        score = methyl_score[patient_id].to_dict()
        if args.layer in PATIENT_SPECIFIC_LAYERS:
            receptors, terminals = Layers.get(args.layer, patient=patient_id)
        else:
            receptors, terminals = Layers.get(args.layer)
        flip = True if args.layer.startswith('terminal') else False
        result = None
        try:
            result = finder.run_average_deregnet(score,
                                                 min_size=10,
                                                 max_size=50,
                                                 time_limit=args.time_limit,
                                                 gap_cut=args.gap_cut,
                                                 abs_values=abs_vals,
                                                 model_sense=minmax,
                                                 receptors=receptors,
                                                 terminals=terminals,
                                                 flip_orientation=flip)
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
