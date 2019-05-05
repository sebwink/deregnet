import os
from biomap import BioMap
from tcga.snv import TcgaSnvData

HGNC = BioMap().get_mapper('hgnc')

DATA_HOME = os.environ.get('DEREGNET_TCGA_DATA_HOME',
                           os.path.join(os.environ['HOME'], 'projects/DeRegNet/deregnet/tcga'))

class Layers:

    @classmethod
    def get(cls, layer, **kwargs):
        if layer == 'vogelstein':
            return cls.vogelstein()
        elif layer == 'vogelstein_oncogenes':
            return cls.vogelstein_oncogenes()
        elif layer == 'vogelstein_tumor_suppressors':
            return cls.vogelstein_tumor_suppressors()
        elif layer == 'vogelstein_as_receptors':
            vogelstein, _ = cls.vogelstein()
            return vogelstein, None
        elif layer == 'vogelstein_as_terminals':
            _, vogelstein = cls.vogelstein()
            return vogelstein, None
        elif layer.startswith('rooted/'):
            genes = layer.split('/')[1]
            return cls.rooted(genes)
        elif layer.startswith('terminal/'):
            genes = layer.split('/')[1]
            return cls.terminal(genes)
        elif layer == 'null':
            return None, None
        elif layer == 'snv':
            snv = cls.snv(**kwargs)
            return snv, snv
        elif layer == 'snv_as_roots':
            snv = cls.snv(**kwargs)
            return snv, None
        elif layer == 'snv_as_terminals':
            snv = cls.snv(**kwargs)
            return snv, None
        elif layer == 'cnv':
            return cls.cnv(**kwargs)
        elif layer == 'genomic':
            return cls.genomic(**kwargs)
        elif layer == 'hla_as_terminals':
            return cls.hla(as_receptors=False, as_terminals=True, **kwargs)
        elif layer == 'svn_to_hla':
            hla, _ = cls.hla(True, False)
            snv, _ = cls.snv(**kwargs)
            return snv, hla
        else:
            raise ValueError

    @classmethod
    def hla(cls, as_receptors, as_terminals):
        hgnc = BioMap().get_mapper('hgnc')
        hla_genes = [ID for ID in hgnc.get_all('symbol') if ID.startswith('HLA')]
        print(hla_genes)
        receptors, terminals = None, None
        if as_receptors:
            receptors = hla_genes
        if as_terminals:
            terminals = hla_genes
        return receptors, terminals


    @classmethod
    def snv(cls, patient, cancer_type):
        snv = TcgaSnvData(cancer_type)
        snvs = snv.non_silent_mutations(patient)
        try:
            snvs = HGNC.map(snvs, FROM='ensembl', TO='entrez')
        except:
            snvs = []
        return [snv for snv in snvs if snv is not None]

    @classmethod
    def cnv(cls, patient):
        return None, None

    @classmethod
    def genomic(cls, patient):
        snvs = cls.snv(patient)
        cnvs = cls.cnv(patient)
        genomic = set(snvs).union(set(cnvs))
        return list(genomic)

    @classmethod
    def rooted(cls, genes):
        return [*HGNC.map(genes.split('_'), FROM='symbol', TO='entrez')], None

    @classmethod
    def terminal(cls, genes):
        return [*HGNC.map(genes.split('_'), FROM='symbol', TO='entrez')], None

    @classmethod
    def vogelstein(cls):
        oncogenes, _ = cls.vogelstein_oncogenes()
        ts, _ = cls.vogelstein_tumor_suppressors()
        vglstn = set(oncogenes).union(set(ts))
        return list(vglstn), list(vglstn)

    @classmethod
    def vogelstein_oncogenes(cls):
        genes = set()
        path = os.path.join(DATA_HOME, 'deregnet_tcga/_layers/vogelstein/oncogenes_vogelstein.txt')
        with open(path, 'r') as fp:
            for line in fp.readlines():
                genes.add(line.strip())
        genes = [ID for ID in HGNC.map(list(genes), FROM='symbol', TO='entrez') if ID]
        return list(genes), list(genes)

    @classmethod
    def vogelstein_tumor_suppressors(cls):
        genes = set()
        path = os.path.join(DATA_HOME, 'deregnet_tcga/_layers/vogelstein/tumor_suppressors_vogelstein.txt')
        with open(path, 'r') as fp:
            for line in fp.readlines():
                genes.add(line.strip())
        genes = [ID for ID in HGNC.map(list(genes), FROM='symbol', TO='entrez') if ID]
        return genes, genes
