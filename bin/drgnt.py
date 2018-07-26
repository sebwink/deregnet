#!/usr/bin/env python3
'''

'''
import argparse

import igraph as ig

from biomap import BioMap

import deregnet.core
import deregnet.script
import deregnet.visual

def define_args(parser):
    deregnet.script.define_shared_args(parser)

    parser.add_argument('--size', default=20, dest='size', type=int,
                        help='Size of the resulting subgraph(s). Default : 20')

    parser.add_argument('--root', default=None, dest='root', type=str,
                        help='Specify a root node. Same id type as graph id type!')

def main():
    # set up argument parser
    parser = argparse.ArgumentParser()
    define_args(parser)
    argparse_args = parser.parse_args()
    # read graph and initialize the subgraph finder
    graph = ig.Graph.Read_GraphML(argparse_args.graph)
    finder = deregnet.SubgraphFinder(graph, argparse_args.graph_id_attr)
    # get the id mapper from BioMap
    id_mapper = BioMap.get_mapper(argparse_args.id_mapper)
    # initialize the argument object ...
    deregnet_args = deregnet.core.AbsoluteDeregnetArguments()
    # ... and populate it
    deregnet.script.populate_shared_args(deregnet_args, argparse_args, id_mapper)
    deregnet_args.size = argparse_args.size
    deregnet_args.root = argparse_args.root
    # find subgraphs
    results = finder.run(deregnet_args)
    # write results
    results.to_graphml(argparse_args.output)
    # TODO: visual reporting

if __name__ == '__main__':
    main()
