from tcga.rnaseq import TcgaRnaSeq
from biomap import BioMap

def prepare_expression_indicator(dataset, threshold, tumor=True):
    id_mapper = BioMap().get_mapper('hgnc')
    if tumor:
        indicator = TcgaRnaSeq(dataset).get_binary_tumor_expression_indicator(threshold=threshold)
    else:
        indicator = TcgaRnaSeq(dataset).get_binary_control_expression_indicator(threshold=threshold)
    indicator.index = [gene.split('.')[0] for gene in indicator.index]
    indicator.index = list(id_mapper.map(list(indicator.index), FROM='ensembl', TO='entrez'))
    entrez = [ID for ID in indicator.index if ID]
    if tumor:
        return indicator.loc[entrez, :]
    else:
        return indicator.loc[entrez]

def get_expression_induced_subgraph(graph, case_id, ei_tumor, ei_control):
    ind_tumor = ei_tumor[case_id].to_dict()
    ind_control = ei_control.to_dict()
    ind = {gene: (indg or ind_control.get(gene, True)) for gene, indg in ind_tumor.items()}
    genes = {g for g, indg in ind.items() if indg}
    return graph.subgraph(graph.vs.select(name_in=genes))
