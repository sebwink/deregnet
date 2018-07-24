import sys
import os
import subprocess

import igraph as ig

from biomap import BioMap
from deregnet.core import SubgraphFinder

sys.path.append('.')

from utils import get_rnaseq_score, get_rnaseq_score_for_patient, get_vogelstein

def main(argv):
    dataset = argv[1]
    graph_path = 'graphs/kegg_hsa_paper.graphml'
    if len(argv) > 2:
        graph_path = argv[2]
    #
    id_mapper = BioMap.get_mapper('hgnc')
    rnaseq_score = get_rnaseq_score(dataset, normalize_wrt='genewise_median')
    rnaseq_score.index = [gene.split('.')[0] for gene in rnaseq_score.index]
    rnaseq_score.index = list(id_mapper.map(list(rnaseq_score.index), FROM='ensembl', TO='entrez'))
    patients = list(rnaseq_score.columns)
    failed = []
    base_path = os.path.join('mode4', dataset)
    vogelstein = id_mapper.map(get_vogelstein(), FROM='symbol', TO='entrez')
    vogelstein = [gene for gene in vogelstein if gene]
    if not os.path.isdir(base_path):
        os.makedirs(base_path)
    for patient_id in patients:
        path = os.path.join(base_path, patient_id)
        nodes_already = set()
        i = 0
        for graphml in os.listdir(path):
            if not graphml.endswith('.graphml'):
                continue
            i += 1
            gf = os.path.join(path, graphml)
            G = ig.Graph.Read_GraphML(gf)
            nodes_already = nodes_already.union({v['name'] for v in G.vs})

        graph = ig.Graph.Read_GraphML(graph_path)
        graph.delete_vertices(graph.vs.select(name_in=nodes_already))
        finder = SubgraphFinder(graph)
        score = get_rnaseq_score_for_patient(patient_id, rnaseq_score)
        try:
            result = finder.run_average_deregnet(score,
                                                 min_size=10,
                                                 gap_cut=0.05,
                                                 time_limit=1200,
                                                 abs_values=True,
                                                 receptors=vogelstein,
                                                 terminals=vogelstein)
        except:
            failed.append(patient_id)
        path = os.path.join(path, 'suboptimal_'+str(i)+'.graphml')
        try:
            result.optimal.write_graphml(path)
        except:
            failed.append(patient_id)
    with open(os.path.join(base_path, 'failed.txt'), 'w') as fp:
        for fail in failed:
            fp.write(fail+'\n')

if __name__ == '__main__':
    main(sys.argv)
