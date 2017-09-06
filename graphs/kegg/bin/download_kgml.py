import os
import sys

import requests

KEGG_REST_ROOT='http://rest.kegg.jp'

def get_pathways(species):
    query = KEGG_REST_ROOT+'/list/pathway/'+species
    response = requests.get(query)
    if response.status_code != 200:
        print('KEGG REST Error.')
        return None
    response = response.content.decode('utf-8')
    pathway_data = [line.strip() for line in response.split('\n') if line]
    return [pathway.split('\t')[0].split(':')[1] for pathway in pathway_data]

def download_kgml(pathway, path):
    query = KEGG_REST_ROOT+'/get/'+pathway+'/kgml'
    response = requests.get(query)
    if response.status_code != 200:
        print('KEGG REST Error.')
        return None
    with open(os.path.join(path, pathway+'.xml'), 'wb') as kgml:
        kgml.write(response.content)

def download_all_kgmls(species, path):
    for pathway in get_pathways(species):
        print('Downloading %s ...' % pathway)
        download_kgml(pathway, path)

if __name__ == '__main__':
    download_all_kgmls(sys.argv[1], sys.argv[2])
