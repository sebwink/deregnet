import os
import pandas as pd

from tcga.snv import TcgaSnvData
from tcga.rnaseq import TcgaRnaSeq
from tcga.methylation import TcgaMethylation450k

def get_somatic_mutation_score(patient_id, mutation_matrix):
    patient_index = mutation_matrix[2].index(patient_id)
    matrix = mutation_matrix[0].tocsc()
    return { gene: matrix[i, patient_index] for i, gene in enumerate(mutation_matrix[1]) }

def get_mutation_matrix(dataset='lihc'):
    data = TcgaSnvData(dataset)
    return data.binary_mutation_matrix_consensus()

__DATAPATH__ = '/home/sewi/projects/DeRegNet/deregnet/data/tcga'

def get_rnaseq_score(dataset='lihc', mode='trinary', compare_to='control'):
    rnaseq = TcgaRnaSeq(dataset)
    if mode=='trinary':
        return rnaseq.get_trinary_log2fold_score(threshold=2, normalize_wrt=compare_to)
    else:
        return rnaseq.get_binary_log2fold_score(threshold=2, normalize_wrt=compare_to)

def get_methylation_score(dataset='lihc', **kwargs):
    methylation = TcgaMethylation(dataset)
    return methylation.get_score(**kwargs)
