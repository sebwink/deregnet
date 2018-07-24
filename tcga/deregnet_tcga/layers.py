import os
from biomap import BioMap

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
        elif layer.startswith('rooted/'):
            genes = layer.split('/')[1]
            return cls.rooted(genes)
        elif layer.startswith('terminal/'):
            genes = layer.split('/')[1]
            return cls.terminal(genes)
        elif layer == 'null':
            return None, None
        elif layer == 'snv':
            return cls.snv(**kwargs)
        elif layer == 'cnv':
            return cls.cnv(**kwargs)
        elif layer == 'genomic':
            return cls.genomic(**kwargs)
        else:
            raise ValueError

    @classmethod
    def snv(cls, patient):
        return None, None

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
