import json
import pandas as pd

class TCGAMetaData(object):
    def __init__(self, metadata_file):
        with open(metadata_file, 'r') as metadata:
             self.metadata = json.load(metadata)

    def search_and_get(self, field, dictionary):
        '''
        Dangerous for fields which appear mutliple times across different nesting levels ...

        Use with caution.
        '''
        if field == 'barcode':
            field = 'entity_submitter_id'
        if field == 'file_id':
            return field, dictionary[field]
        for key, value in dictionary.items():
            k, v = None, None
            if key == field and not isinstance(value, dict) and not isinstance(value,list):
                return field, value
            elif isinstance(value, dict):
                 k, v = self.search_and_get(field, value)
            elif isinstance(value, list):
                 for elmt in value:
                     k, v = self.search_and_get(field, elmt)
                     if k is not None: break
            if k is not None:
                return field, v

        return None, None

    def search(self, field, dictionary):
        k, v = self.search_and_get(field, dictionary)
        if k is None:
            return False
        return True

    def create_table(self, columns, renamer = {}, transformer = {}):
        columns = set(columns).union({'file_id'})
        for col in columns:
            if col not in renamer.keys():
                renamer[col] = col
            if col not in transformer.keys():
                transformer[col] = lambda c: c

        table_dict = {}
        for field in columns:
            table_dict[renamer[field]] = [transformer[field](self.search_and_get(field, file)[1]) for file in self.metadata]
        table = pd.DataFrame(table_dict)
        # table.set_index('file_id', inplace = True)
        return table
