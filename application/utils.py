import os
import pandas as pd
from tcga.analysis.wxs_somatic_masked import TcgaWxsMaskedMutationData

def get_somatic_mutation_score(patient_id, mutation_matrix):
    patient_index = mutation_matrix[2].index(patient_id)
    matrix = mutation_matrix[0].tocsc()
    return { gene: matrix[i, patient_index] for i, gene in enumerate(mutation_matrix[1]) }

def get_mutation_matrix(dataset='lihc'):
    data = TcgaWxsMaskedMutationData(dataset)
    return data.binary_mutation_matrix_consensus()

__DATAPATH__ = '/home/sebastian/wrk/deregnet/deregnet/data/tcga'

def get_rnaseq_score(dataset='lihc'):
    score_path = os.path.join(__DATAPATH__, dataset, 'rnaseq/trinary_log2fold_score.csv')
    return pd.read_csv(score_path)

def get_rnaseq_score_for_patient(patient_id, score_table):
    return score_table[patient_id].to_dict()

def get_vogelstein(path='/home/sebastian/wrk/deregnet/deregnet/data/cancer_genes_vogelstein.txt'):
    genes = set()
    with open(path, 'r') as fp:
        for line in fp.readlines():
            genes.add(line.strip())
    return list(genes)
