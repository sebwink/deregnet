#!/share/opt/x86_64_sl7/anaconda-python-3.6/bin/python

import os
import subprocess
import time
import argparse
import igraph

import deregnet

from biograph.mapping.gene import HGNCMapper

def define_args(parser):
    parser.add_argument('--absolute-values', default = False, dest = 'abs', type = bool,
                        help = 'Whether to take absolute values of the scores. Default: False, use like ... --absolute-values True ...')
    parser.add_argument('--graph', dest = 'graph', type = str, required=True,
                        help = 'A graphml file containing the graph you want to run drgnt with.')
    parser.add_argument('--scores', metavar = 'SCORE_FILE', dest = 'scores', type = str, required=True,
                        help = 'A file containing the scores. See further options below.')
    parser.add_argument('--score-column', default = 'score', dest = 'score_col', type = str, 
                        help = 'Column name of (gene) id in your score file. Default: score')
    parser.add_argument('--id-column', default = 'id', dest = 'id_col', type = str,
                        help = 'Column name of (gene) id in your score file. Default: id')
    parser.add_argument('--sep', default = 'tab', dest = 'sep', type = str,
                        help = 'The column seperator in your score file. Options: comma, tab. Default: tab')
    parser.add_argument('--score-id-type', default = None, dest = 'score_id_type', type = str,
                        help = 'Which id type do you have in your score file? Options: entrez, uniprot, ensembl, symbol. Default: same as graph id type')
    parser.add_argument('--graph-id-type', default = 'entrez', dest = 'graph_id_type', type = str,
                        help = 'Which id type does the graph have? Options: entrez, uniprot, ensembl, symbol. Default:entrez')
    parser.add_argument('--graph-id-attr', default = 'name', dest = 'graph_id_attr', type = str,
                        help = 'Node attribute which contains the relevant id in the graphml. Default: name')
    parser.add_argument('--species', default = 'hsa', dest = 'species', type = str,
                        help = 'Species your id\'s refer to. KEGG code: hsa, eco, sce, etc... Default: hsa (Homo Sapiens)')
    parser.add_argument('--min-size', default = None, dest = 'min_size', type = str,
                        help = 'Minimal size of the resulting subgraph(s). Default : 20')
    parser.add_argument('--max-size', default = None, dest = 'max_size', type = str,
                        help = 'Maximal size of the resulting subgraph(s). Default : 50')
    parser.add_argument('--root', default = None, dest = 'root', type = str,
                        help = 'Specify a root node. Same id type as graph id type!')
    parser.add_argument('--suboptimal', default = None, dest = 'suboptimal', type = str,
                        help = 'Number of suboptimal subgraphs you want to find. (Increases runtime)')
    parser.add_argument('--max-overlap-percentage', default = None, dest = 'max_overlap', type = str,
                        help = 'How much can suboptimal subgraphs overlap with already found subgraphs. Default: 0')
    parser.add_argument('--gap-cut', default = None, dest = 'gap', type = str,
                        help = 'Stop optimization prematurely if current solution within GAP of optimal solution. Default: None')
    parser.add_argument('--time-limit', default = None, dest = 'time_limit', type = str,
                        help = 'Set a time limit in seconds. Default: None')
    parser.add_argument('--model_sense', default = None, dest = 'model_sense', type = str,
                        help = 'Model sense. \noptions: max, min. Default: max')
    parser.add_argument('--output-path', default = os.getcwd(), dest = 'output', type = str,
                        help = 'Folder to which output is written. (Does not have to exist.) Default : cwd')
    parser.add_argument('--flip-orientation', default = False, dest='flip', type = bool,
                        help = 'Set --flip-orientation True when you want to flip the orientation. Default: False')
    parser.add_argument('--algorithm', default = 'dta', dest = 'algorithm', type = str, 
                        help = 'Algorithm to use: dta, gcc, ovt.')
    parser.add_argument('--receptors-file', default=None, dest='receptors', type=str,
                        help = 'Specify path of file containing the receptor nodes')
    parser.add_argument('--terminals-file', default=None, dest='terminals', type=str,
                        help = 'Specify path of file containing the terminal nodes')
    

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
                                         HGNCMapper)

    drgnt_arg_dict = {
                    '--graph' : tmp_files.graph,
                    '--score' : tmp_files.scores,
                    '--output-dir' : os.path.join(tmp_files.path, 'subgraphs'),
                    '--gap-cut' : args.gap,
                    '--min-size' : args.min_size,
                    '--max-size': args.max_size,
                    '--suboptimal' : args.suboptimal,
                    '--max-overlap-percentage' : args.max_overlap,
                    '--root' : args.root,
                    '--time-limit' : args.time_limit,
                    '--model-sense' : args.model_sense,
                    '--algorithm' : args.algorithm,
                    '--receptors-file' : args.receptors,
                    '--terminals-file' : args.terminals
                 }

    drgnt_args = []
    for arg in drgnt_arg_dict:
        if drgnt_arg_dict[arg] is not None:
             drgnt_args += [arg, drgnt_arg_dict[arg]]

    if args.flip:
        drgnt_args += ['--flip-orientation']
    if args.abs:
        drgnt_args += ['--absolute-values']

    subprocess.call(['avgdrgnt'] + drgnt_args)

    if not os.path.isdir(args.output):
        os.makedirs(args.output)
    deregnet.output2graphml(graph, tmp_files, args.output)

if __name__ == '__main__':
    main()
