import os
import glob

import numpy as np
import scipy.sparse
import pandas as pd

__mypath__ = os.path.dirname(os.path.abspath(__file__))
__tcga_home__ = __mypath__

class TcgaWxsMaskedMutationData(object):

    def __init__(self, cancer_type = 'lihc', path = None):
        if path is None:
            path = os.path.join(__tcga_home__, cancer_type + '/wxs_somatic_masked')
        print(path)
        self.path = path
        maf_files = []

        for directory in os.listdir(path):
            path2dir = os.path.join(path, directory)
            if os.path.isdir(path2dir):
                maf_files.append(glob.glob(os.path.abspath(os.path.join(path2dir, '*.somatic.maf.gz')))[0])
        methods = ['muse', 'varscan', 'mutect', 'somaticsniper']
        self.datasets = {
                          method : {'maf': maf_file, 'data': None }

                          for method in methods
                          for maf_file in maf_files
                          if method in maf_file 
                         }

    def load(self, method = 'varscan'):
        if self.datasets[method]['data'] is None:
            self.datasets[method]['data'] = pd.read_table(self.datasets[method]['maf'],
                                                          compression = 'gzip',
                                                          comment = '#',
                                                          low_memory = False)

    def delete(self, method = 'varscan'):
        self.datasets[method]['data'] = None

    def column_filter(self, method = 'varscan',
                      exclude = {},
                      **kwargs):

        self.load(method)
        df = self.datasets[method]['data']
        if not kwargs:
            return df
        column_filters = []
        for column, values in kwargs.items():
            exclude_values = exclude.get(column, True)
            if exclude_values:
                column_filter = (df[column] != values[0])
            else:
                column_filter = (df[column] == values[0])
            for i in range(1, len(values)):
                if exclude_values:
                    column_filter = column_filter | (df[column] != values[i])
                else:
                    column_filter = column_filter | (df[column] == values[i])
            column_filters.append(column_filter)
        filter_expr = column_filters[0]
        for i in range(1, len(column_filters)):
            filter_expr = filter_expr & column_filters[i]
        return df[filter_expr]

    def get_mutations(self, method = 'varscan',
                      exclude = {},
                      **kwargs):
        self.load(method)
        filtered_data = self.column_filter(method, exclude, **kwargs)
        filtered_data = filtered_data.dropna(0, subset = ['Gene'])
        return {
                 (t.Gene, t.Tumor_Sample_Barcode.split('-')[2])
                  for t in filtered_data[['Gene', 'Tumor_Sample_Barcode']].itertuples()
               }

    def get_matrix_from_mutations(self, mutations):
        genes = list({mutation[0] for mutation in mutations})
        patients = list({mutation[1] for mutation in mutations})
        mutation_matrix = scipy.sparse.dok_matrix((len(genes), len(patients)), dtype = np.int8)
        for mutation in mutations:
            mutation_matrix[genes.index(mutation[0]), patients.index(mutation[1])] = 1
        return scipy.sparse.coo_matrix(mutation_matrix), genes, patients

    def binary_mutation_matrix(self, method = 'varscan',
                               exclude = {},
                               **kwargs):
        mutations = self.get_mutations(method, exclude, **kwargs)
        return self.get_matrix_from_mutations(mutations)

    def binary_mutation_matrix_consensus(self, methods = ['varscan', 'muse', 'mutect', 'somaticsniper'],
                                         consensus_policy = 'intersection',
                                         exclude = {},
                                         **kwargs):
        method2mutations = { method : self.get_mutations(method, exclude, **kwargs) for method in methods }
        if isinstance(consensus_policy, str):
            consensus_policy = eval('self.' + consensus_policy + '_consensus')
        mutations = consensus_policy(method2mutations, methods)
        return self.get_matrix_from_mutations(mutations)

    def intersection_consensus(self, method2mutations, methods):
        mutations = method2mutations[methods[0]]
        for method in method2mutations:
            if method != methods[0]:
                mutations = mutations.intersection(method2mutations[method])
        return mutations

    def union_consensus(self, method2mutations, methods):
        mutations = set()
        for method in method2mutations:
            mutations = mutations.union(method2mutations[method])
        return mutations
