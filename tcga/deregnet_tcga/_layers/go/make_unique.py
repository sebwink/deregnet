import sys
from biomap import BioMap

hgnc = BioMap().get_mapper('hgnc')

FILE = sys.argv[1]
with open(FILE, 'r') as fp:
    uniprot = fp.read().split('\n')
    uniprot = [ID for ID in uniprot if ID]
    uniprot = list(set(uniprot))
    uniprot = [ID.split(':')[1] for ID in uniprot]

symbol = hgnc.map(uniprot, FROM='uniprot')
symbol = [s for s in symbol if s]

with open(FILE, 'w') as fp:
    fp.write(symbol[0])
    for s in symbol[1:]:
        fp.write('\n'+s)
