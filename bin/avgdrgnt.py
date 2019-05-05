#!/usr/bin/env python3
'''

'''
import argparse

import igraph as ig

# from biomap import BioMap

import deregnet.core
import deregnet.script
import deregnet.visual

def define_args(parser):
    deregnet.script.define_shared_args(parser)

    parser.add_argument('--min-size', default=15, dest='min_size', type=int,
                        help='Minimal size of the resulting subgraph(s). Default : 15')

    parser.add_argument('--max-size', default=15, dest='max_size', type=int,
                        help='Maximal size of the resulting subgraph(s). Default : 15')

    parser.add_argument('--min-num-terminals', default=0, dest='min_num_terminals', type=int,
                        help='Minimum number of terminals in the resulting subgraph(s). Default : 0')

    parser.add_argument('--algorithm', default='gcc', dest='algorithm', type=str,
                        choices=['GeneralizedCharnesCooper', 'Dinkelbach', 'ObjectiveVariableTransform'],
                        help='Algorithm to use to solve the fractional integer programming problem.'\
                             'Default: GeneralizedCharnesCooper.')

    deregnet.script.define_geneset_arg(parser, 'receptor')
    deregnet.script.define_geneset_arg(parser, 'terminal')

def main():
    # set up argument parser
    parser = argparse.ArgumentParser()
    define_args(parser)
    argparse_args = parser.parse_args()
    # read graph and initialize the subgraph finder
    graph = ig.Graph.Read_GraphML(argparse_args.graph)
    finder = deregnet.SubgraphFinder(graph, argparse_args.graph_id_attr)
    # get the id mapper from BioMap
    # id_mapper = BioMap().get_mapper(argparse_args.id_mapper)
    id_mapper = None
    # initialize the argument object ...
    deregnet_args = deregnet.core.AverageDeregnetArguments()
    # ... and populate it
    deregnet.script.populate_shared_args(deregnet_args, argparse_args, id_mapper)
    deregnet_args.min_size = argparse_args.min_size
    deregnet_args.max_size = argparse_args.max_size
    deregnet_args.min_num_terminals = argparse_args.min_num_terminals
    deregnet_args.algorithm = argparse_args.algorithm
    deregnet_args.receptors = deregnet.script.parse_geneset('receptor', argparse_args, id_mapper)
    deregnet_args.terminals = deregnet.script.parse_geneset('terminal', argparse_args, id_mapper)
    # find subgraphs
    results = finder.run(deregnet_args)
    # write results
    results.to_graphml(argparse_args.output)
    # TODO: visual reporting

if __name__ == '__main__':
    main()
