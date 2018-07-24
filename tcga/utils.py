import os
import pandas as pd

from tcga.snv import TcgaSnvData
from tcga.rnaseq import TcgaRnaSeq

def get_somatic_mutation_score(patient_id, mutation_matrix):
    patient_index = mutation_matrix[2].index(patient_id)
    matrix = mutation_matrix[0].tocsc()
    return { gene: matrix[i, patient_index] for i, gene in enumerate(mutation_matrix[1]) }

def get_mutation_matrix(dataset='lihc'):
    data = TcgaSnvData(dataset)
    return data.binary_mutation_matrix_consensus()

__DATAPATH__ = '/home/sewi/projects/DeRegNet/deregnet/data/tcga'

def get_rnaseq_score(dataset='lihc', mode='trinary', normalize_wrt='control'):
    rnaseq = TcgaRnaSeq(dataset)
    if mode=='trinary':
        return rnaseq.get_trinary_log2fold_score(threshold=2, normalize_wrt=normalize_wrt)
    else:
        return rnaseq.get_binary_log2fold_score(threshold=2, normalize_wrt=normalize_wrt)

def get_rnaseq_score_for_patient(patient_id, score_table):
    return score_table[patient_id].to_dict()

def get_vogelstein(path='/home/sewi/projects/DeRegNet/deregnet/data/cancer_genes_vogelstein.txt'):
    genes = set()
    with open(path, 'r') as fp:
        for line in fp.readlines():
            genes.add(line.strip())
    return list(genes)
