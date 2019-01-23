import os
import pandas as pd

from biomap import BioMap

from tcga.snv import TcgaSnvData
from tcga.cnv import CnvData
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

def get_global_rnaseq_score(dataset, score):
    rnaseq = TcgaRnaSeq(dataset)
    hgnc = BioMap().get_mapper('hgnc')
    de = rnaseq.get_deseq2_differential_expression()
    de['entrez'] = hgnc.map([ID.split('.')[0] for ID in de.index], FROM='ensembl', TO='entrez')
    if score == 'log2fold':
        scores = {de.ix[i,'entrez']: de.ix[i,'log2FoldChange'] for i in range(len(de)) if de.ix[i,'entrez']}
    else:
        pass
    return scores

def get_cnv_score(dataset='lihc', score='tumor_segment_mean'):
    cnv = CnvData(dataset)
    return cnv.genelevel_score(score)

def get_methylation_score(dataset='lihc', **kwargs):
    methylation = TcgaMethylation450k(dataset)
    return methylation.get_score(**kwargs)

def get_rnaseq_data(ctype='lihc'):
    rnaseq = TcgaRnaSeq('lihc')
    log2fold_wrt_control = rnaseq.get_log2folds(normalize_wrt='control')
    counts = rnaseq.get_deseq2_normalized_counts(patient_barcode=True)
    counts = counts.median(axis=1, level=0)
    log2fold_wrt_control = log2fold_wrt_control.median(axis=1, level=0)
    counts.index = [ID.split('.')[0] for ID in counts.index]
    log2fold_wrt_control.index = [ID.split('.')[0] for ID in log2fold_wrt_control.index]
    return counts, log2fold_wrt_control
