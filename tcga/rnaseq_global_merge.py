import os
import sys
import glob
import igraph as ig

from tcga.rnaseq import TcgaRnaSeq
from biomap import BioMap

rnaseq = TcgaRnaSeq(sys.argv[1])
de = rnaseq.get_deseq2_differential_expression()
hgnc = BioMap().get_mapper('hgnc')
de['symbol'] = hgnc.map([ID.split('.')[0] for ID in de.index], FROM='ensembl')


KEGG = ig.Graph.Read_GraphML('graph/kegg_hsa.graphml')
root = os.path.join('.', 'rnaseq', 'global')
ROOT = os.path.join(root, sys.argv[1], sys.argv[2])

mode = sys.argv[3] if len(sys.argv) == 4 else None

nodes = set()

for root, dirs, files in os.walk(ROOT):
    if mode is not None and mode not in root:
        continue
    graphmls = glob.glob(os.path.join(root, '*.graphml'))
    for graphml in graphmls:
        subgraph = ig.Graph.Read_GraphML(graphml)
        nodes = nodes.union(set(subgraph.vs['symbol']))

graph = KEGG.subgraph(KEGG.vs.select(symbol_in=nodes))
de = de[de['symbol'].isin(nodes)]
for v in graph.vs:
    row = de[de['symbol'] == v['symbol']]
    if len(row) > 0:
        for col in de.columns:
            v[col] = row.ix[0, col]
        if sys.argv[2] == 'log2fold':
            v['deregnet_score'] = v['log2FoldChange']

if mode is None:
    graph.write_graphml('./rnaseq/global/'+sys.argv[1]+'/'+sys.argv[2]+'/merge.graphml')
else:
    graph.write_graphml('./rnaseq/global/'+sys.argv[1]+'/'+sys.argv[2]+'/merge_'+mode+'.graphml')
