import sys
import os

import igraph as ig

from biomap import BioMap
from deregnet.core import SubgraphFinder
from utils import *

def main(argv):
    dataset = argv[1]
    graph_path = 'graphs/kegg_hsa_paper.graphml'
    if len(argv) > 2:
        graph_path = argv[2]
    graph = ig.Graph.Read_GraphML(graph_path)
    finder = SubgraphFinder(graph)
    #
    id_mapper = BioMap.get_mapper('hgnc')
    mutation_data = list(get_mutation_matrix(dataset))
    mutation_data[1] = id_mapper.map(mutation_data[1], FROM='ensembl', TO='entrez')
    patients = mutation_data[2]
    failed = []
    base_path = os.path.join('mode2', dataset)
    os.makedirs(base_path)
    for patient_id in patients:
        score = get_somatic_mutation_score(patient_id, mutation_data)
        result = finder.run_average_deregnet(score, min_size=10, time_limit=1200)
        path = os.path.join(base_path, patient_id)
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
