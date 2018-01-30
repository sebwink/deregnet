from tcga.analysis.wxs_somatic_masked import TcgaWxsMaskedMutationData

def get_somatic_mutation_score(patient_id, mutation_matrix):
    patient_index = mutation_matrix[2].index(patient_id)
    matrix = mutation_matrix[0].tocsc()
    return { gene: matrix[i, patient_index] for i, gene in enumerate(mutation_matrix[1]) }

def get_mutation_matrix(dataset='lihc'):
    data = TcgaWxsMaskedMutationData(dataset)
    return data.binary_mutation_matrix_consensus()
