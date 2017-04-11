#!/opt/anaconda/3.5/bin/python

import os
import subprocess
import time
import argparse
import igraph

from biograph.mapping.gene import GeneIdMapper

import deregnet

def define_args(parser):
    parser.add_argument('--graph', dest = 'graph', type = str)
    parser.add_argument('--scores', dest = 'scores', type = str)
    parser.add_argument('--gap-cut', default = None, dest = 'gap', type = str)
    parser.add_argument('--score-column', default = 'score', dest = 'score_col', type = str)
    parser.add_argument('--id-column', default = 'id', dest = 'id_col', type = str)
    parser.add_argument('--sep', default = ',', dest = 'sep', type = str)
    parser.add_argument('--score-id-type', default = None, dest = 'score_id_type', type = str)
    parser.add_argument('--graph-id-type', default = 'entrez', dest = 'graph_id_type', type = str)
    parser.add_argument('--graph-id-attr', default = 'name', dest = 'graph_id_attr', type = str)
    parser.add_argument('--species', default = 'hsa', dest = 'species', type = str)
    parser.add_argument('--size', default = None, dest = 'size', type = str)
    parser.add_argument('--suboptimal', default = None, dest = 'suboptimal', type = str)
    parser.add_argument('--max-overlap-percentage', default = None, dest = 'max_overlap', type = str)
    parser.add_argument('--root', default = None, dest = 'root', type = str)
    parser.add_argument('--time-limit', default = None, dest = 'time_limit', type = str)
    parser.add_argument('--model_sense', default = None, dest = 'model_sense', type = str)
    parser.add_argument('--output-path', default = os.getcwd(), dest = 'output', type = str)

def parse_sep(sep):
    if sep == 'tab':
        return '\t'
    return ','

def main():
    parser = argparse.ArgumentParser()
    define_args(parser)
    args = parser.parse_args()

    graph = igraph.Graph.Read(args.graph, 'graphml')
    if args.score_id_type is None:
        args.score_id_type = args.graph_id_type

    tmp_files = deregnet.write_tmp_files(graph,
                                         args.graph_id_attr,
                                         args.graph_id_type,
                                         args.scores,
                                         args.score_id_type,
                                         parse_sep(args.sep),
                                         args.id_col,
                                         args.score_col,
                                         args.species,
                                         GeneIdMapper)

    drgnt_arg_dict = {
                    '--graph' : tmp_files.graph,
                    '--score' : tmp_files.scores,
                    '--output-dir' : os.path.join(tmp_files.path, 'subgraphs'),
                    '--gap-cut' : args.gap,
                    '--size' : args.size,
                    '--suboptimal' : args.suboptimal,
                    '--max-overlap-percentage' : args.max_overlap,
                    '--root' : args.root,
                    '--time-limit' : args.time_limit,
                    '--model-sense' : args.model_sense
                 }

    drgnt_args = []
    for arg in drgnt_arg_dict:
        if drgnt_arg_dict[arg] is not None:
            drgnt_args += [arg, drgnt_arg_dict[arg]]

    subprocess.call(['drgnt'] + drgnt_args)

    if not os.path.isdir(args.output):
        os.makedirs(args.output)
    deregnet.output2graphml(graph, tmp_files, args.output)

if __name__ == '__main__':
    main()
