#!/usr/bin/env python3

import os
import sys
import argparse
import subprocess

from biomap import BioMap

from deregnet_tcga.layers import Layers

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
    parser.add_argument('-s', '--script', type=str, default=None,
                        help='Which script to run.')
    return parser.parse_args()

args = parse_args()

_, terminals = Layers.get(args.layer)

print(terminals)

terminals_entrez = BioMap().get_mapper('hgnc').map(terminals, FROM='symbol', TO='entrez')
terminals_entrez = [ID for ID in terminals_entrez if ID]

for gene in terminals:
    gene_entrez = BioMap().get_mapper('hgnc').map(gene, TO='entrez')
    print(gene)
    print(args.script)
    run = ['python3',
           args.script,
           '--layer',
           'terminal/'+gene,
           '-g',
           str(args.gap_cut),
           '-t',
           str(args.time_limit),
           '-d',
           args.dataset,
           '-m',
           args.mode,
           '--min-size',
           '5', 
           '--exclude',
           ','.join([ID for ID in set(terminals_entrez).difference(gene_entrez)])]
    subprocess.run(run)
