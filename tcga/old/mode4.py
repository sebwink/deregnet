import sys
import os

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
    graph = ig.Graph.Read_GraphML(graph_path)
    finder = SubgraphFinder(graph)
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
        if os.path.isfile(os.path.join(path, 'optimal.graphml')):
            continue
        score = get_rnaseq_score_for_patient(patient_id, rnaseq_score)
        print(score)
        try:
            result = finder.run_average_deregnet(score,
                                                 min_size=10,
                                                 time_limit=1200,
                                                 abs_values=True,
                                                 receptors=vogelstein,
                                                 terminals=vogelstein)
        except:
            failed.append(patient_id)
        if not os.path.isdir(path):
            os.makedirs(path)
        try:
            result.to_graphml(path)
        except:
            failed.append(patient_id)
    with open(os.path.join(base_path, 'failed.txt'), 'w') as fp:
        for fail in failed:
            fp.write(fail+'\n')

if __name__ == '__main__':
    main(sys.argv)
