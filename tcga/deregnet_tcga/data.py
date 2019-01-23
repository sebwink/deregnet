import os
import pandas as pd

DATA_PATH = os.path.join(os.path.dirname(__file__), '../data/lihc/clinical')

def lihc_clinical_data(path=DATA_PATH):
    clinical_data = os.path.join(path, 'clinical.tsv')
    exposure_data = os.path.join(path, 'exposure.tsv')
    return pd.read_table(clinical_data), pd.read_table(exposure_data)
