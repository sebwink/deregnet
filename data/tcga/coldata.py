import sys
import tcga.analysis.metadata as metadata

def sample_type_id_summary(sample_type_id):
    if sample_type_id in [''.join(['0', str(i)]) for i in range(1,10)]:
        return 'tumor'
    elif sample_type_id in [''.join(['1', str(i)]) for i in range(10)]:
        return 'normal'
    return 'special'

def construct(path2metadata):
    md = metadata.TCGAMetaData(path2metadata)
    return md.create_table(['sample_type_id'], {'sample_type_id':'condition'}, {'sample_type_id': sample_type_id_summary})

if __name__ == '__main__':
    coldata = construct(sys.argv[1])
    coldata.to_csv('coldata.tsv.gz', sep='\t', compression='gzip', index=False)

