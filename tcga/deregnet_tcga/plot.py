from tcga.snv import TcgaSnvData
from biomap import BioMap
import seaborn as sns

def patient_to_genes_affected_by_snv(cancer_type, id_type=None, **kwargs):
    snv_data = TcgaSnvData(cancer_type)
    mutations = snv_data.get_mutations(**kwargs)
    patients = {mutation[1] for mutation in mutations}
    mapping = {
               patient: {mutation[0] for mutation in mutations if mutation[1] == patient}
               for patient in patients
              }
    if id_type is not None:
        hgnc = BioMap().get_mapper('hgnc')
        mapping = {patient: set(hgnc.map(list(genes), FROM='ensembl', TO=id_type))
                   for patient, genes in mapping.items()}
        mapping = {patient: {gene for gene in genes if gene}
                   for patient, genes in mapping.items()}
    return mapping

def get_mutation_color(matrix, patient2genes, genes, color='red'):
    genes = genes if isinstance(genes, list) else [genes]
    colors = []
    color_palette = [color]
    for patient in matrix.columns:
        color_set = False
        for i, gene in enumerate(genes):
            if gene in patient2genes.get(patient, {}):
                colors.append(color_palette[i])
                color_set = True
                continue
        if not color_set:
            colors.append('grey')
    return colors

def make_mutation_clustermap(distdata, patient2genes, gene1, gene2=None, color1='yellow', color2='blue',vmax=0.25, **kwargs):
    colors1 = get_mutation_color(distdata, patient2genes, gene1, color1)
    if gene2:
        colors2 = get_mutation_color(distdata, patient2genes, gene2, color2)
    sns.clustermap(distdata,
                   row_colors=colors1,
                   col_colors=colors2 if gene2 else colors1,
                   xticklabels=False,
                   yticklabels=False,
                   vmax=vmax,
                   **kwargs)
