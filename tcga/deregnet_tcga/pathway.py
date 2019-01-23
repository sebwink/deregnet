import pandas as pd
from biomap import BioMap
from datautils.pathview import PathViewR

def make_contrast(cases, log2fold):
    contrast = pd.DataFrame()
    contrast['case'] = log2fold[cases].mean(axis=1)
    not_cases = list(set(log2fold.columns) - set(cases))
    contrast['not_case'] = log2fold[not_cases].mean(axis=1)
    contrast.dropna()
    hgnc = BioMap().get_mapper('hgnc')
    contrast['entrez'] = hgnc.map(list(contrast.index), FROM='ensembl', TO='entrez')
    contrast = contrast.groupby('entrez').mean()
    return contrast

def pathview(subgraphs, log2fold_wrt_control, genes, pathway, path, description):
    cases = subgraphs.with_any_of(genes)
    contrast = make_contrast(cases, log2fold_wrt_control)
    return PathViewR.run(contrast, pathway, path, description)

def pathview_cases(cases, log2fold_wrt_control, pathway, path, description):
    contrast = make_contrast(cases, log2fold_wrt_control)
    return PathViewR.run(contrast, pathway, path, description)