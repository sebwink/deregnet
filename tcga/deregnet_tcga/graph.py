from tcga.rnaseq import TcgaRnaSeq
from biomap import BioMap

def prepare_expression_indicator(dataset, threshold):
    id_mapper = BioMap().get_mapper('hgnc')
    indicator = TcgaRnaSeq(dataset).get_binary_tumor_expression_indicator(threshold=threshold) 
    indicator.index = [gene.split('.')[0] for gene in indicator.index]
    indicator.index = list(id_mapper.map(list(indicator.index), FROM='ensembl', TO='entrez'))
    entrez = [ID for ID in indicator.index if ID]
    return indicator.loc[entrez, :]

def get_expression_induced_subgraph(graph, case_id, ei):
    ind = ei[case_id].to_dict()
    genes = {g for g, i in ind.items() if i == 1}
    return graph.subgraph(graph.vs.select(name_in=genes))

