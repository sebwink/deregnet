'''

'''
import os
import math

import pandas as pd

import biomap.utils.geneset as geneset

def define_geneset_arg(parser, layer_name):
    # --<layer_name>-file (GMT or GRP format)
    parser.add_argument('--'+layer_name+'-file', dest=layer_name+'_file', type=str, action='append', default=None,
                        help='Path to GMT or GRP file containing genes defining the '+layer_name+' layer.')
    # --<layer_name>-genesets (comma seperated list of geneset names, only applicable if GMT file provided)
    parser.add_argument('--'+layer_name+'-genesets', dest=layer_name+'_genesets', type=str, action='append', default=None,
                        help='Comma seperated list of geneset names for '+layer_name+' layer,'\
                             'only applicable if GMT file provided.')
    # --<layer_name> (comma seperated list)
    parser.add_argument('--'+layer_name, dest=layer_name, type=str, default=None,
                        help='Comma seperated list of IDs defining the '+layer_name+' layer.')
    # --<layer_name>-id-type
    parser.add_argument('--'+layer_name+'-id-type', dest=layer_name+'_id_type', type=str, default=None,
                        help='Id-type for '+layer_name+' layer genesets. Options: all supported by chosen biomap mapper')

def define_shared_args(parser):

    define_geneset_arg(parser, 'include')
    define_geneset_arg(parser, 'exclude')

    parser.add_argument('--debug', action='store_true', help='Debug underlying C++ code with gdb.')

    parser.add_argument('--absolute-values', action='store_true', dest='abs_values',
                        help='Whether to take absolute values of the scores.')

    parser.add_argument('--graph', dest='graph', type=str, required=True,
                        help='A graphml file containing the graph you want to run drgnt with.')

    parser.add_argument('--scores', metavar='SCORE_FILE', dest='scores', type=str, required=True,
                        help='A text file containing the scores. See further options below.')

    parser.add_argument('--default-score', default=0.0, dest='default_score', type=float,
                        help='The score of nodes in the graph which are not scored in your score file. Default: 0.0')

    parser.add_argument('--score-column', default='score', dest='score_col', type=str,
                        help='Column name of (gene) id in your score file. Default: score')

    parser.add_argument('--score-file-without-header', action='store_false', dest='has_header',
                        help='Flag to indicate whether the score file has a header or not.')

    parser.add_argument('--id-column', default='id', dest='id_col', type=str,
                        help='Column name of (gene) id in your score file. Default: id')

    parser.add_argument('--sep', default='\t', dest='sep', type=str,
                        help='The column seperator in your score file.'\
                              'Options: comma, tab. Default: \\t')

    parser.add_argument('--biomap-mapper', default='hgnc', dest='id_mapper', type=str,
                        help='biomap mapper you want to use for id mapping. Default: hgnc')

    parser.add_argument('--score-id-type', default=None, dest='score_id_type', type=str,
                        help='Which id type do you have in your score file? Options: all those'\
                             'supported by the biomap mapper you chose or unspecified. Default: same as graph id type')

    parser.add_argument('--graph-id-type', default=None, dest='graph_id_type', type=str,
                        help='Which id type does the graph have? Options: all those supported'\
                              'by the biomap mapper you chose or unspecified. Default: unspecifed i.e. None')

    parser.add_argument('--graph-id-attr', default='name', dest='graph_id_attr', type=str,
                        help = 'Node attribute which contains the relevant id in the graphml. Default: name')

    parser.add_argument('--suboptimal', default=0, dest='suboptimal', type=int,
                        help = 'Number of suboptimal subgraphs you want to find. (Increases runtime)')

    parser.add_argument('--max-overlap-percentage', default=0, dest='max_overlap', type=float,
                        help = 'How much can suboptimal subgraphs overlap with already found subgraphs. Default: 0')

    parser.add_argument('--gap-cut', default=None, dest='gap_cut', type=float,
                        help = 'Stop optimization prematurely if current solution within GAP of optimal solution. Default: None')

    parser.add_argument('--time-limit', default=None, dest='time_limit', type=int,
                        help = 'Set a time limit in seconds. Default: None')

    parser.add_argument('--model_sense', default='max', dest='model_sense', type=str, choices=['min', 'max'],
                        help='Model sense. Default: max')

    parser.add_argument('--output-path', default=os.getcwd(), dest='output', type=str,
                        help='Folder to which output is written. (Does not have to exist.) Default : cwd')

    parser.add_argument('--flip-orientation', action='store_true',
                        help = 'Set --flip-orientation when you want to flip the orientation of the underlying graph.')

def populate_shared_args(deregnet_args, argparse_args, id_mapper):
    '''
    '''
    deregnet_args.scores = parse_scores(argparse_args.scores,
                                        argparse_args.score_col,
                                        argparse_args.id_col,
                                        argparse_args.sep,
                                        argparse_args.has_header,
                                        id_mapper,
                                        argparse_args.score_id_type,
                                        argparse_args.graph_id_type)
    deregnet_args.default_score = argparse_args.default_score
    deregnet_args.included_nodes = parse_geneset('include', argparse_args, id_mapper)
    deregnet_args.excluded_nodes = parse_geneset('exclude', argparse_args, id_mapper)
    deregnet_args.flip_orientation = argparse_args.flip_orientation
    deregnet_args.model_sense = argparse_args.model_sense
    deregnet_args.time_limit = argparse_args.time_limit
    deregnet_args.gap_cut = argparse_args.gap_cut
    deregnet_args.max_overlap = argparse_args.max_overlap
    deregnet_args.suboptimal = argparse_args.suboptimal
    deregnet_args.debug = argparse_args.debug
    deregnet_args.abs_values = argparse_args.abs_values

def parse_scores(score_file,
                 score_column,
                 id_column,
                 col_sep='\t',
                 has_header=True,
                 id_mapper=None,
                 score_id_type=None,
                 graph_id_type=None):
    '''

    '''
    header = {'header':None} if not has_header else {}
    if not has_header:
        score_column, id_column = int(score_column), int(id_column)
    if isinstance(score_file, pd.DataFrame):
        df = score_file
    else:
        df = pd.read_table(score_file, sep=col_sep, dtype={'id_column': str}, **header)
    if score_id_type != graph_id_type:
        ids = list(df[id_column])
        ids = id_mapper.map(ids, score_id_type, graph_id_type)
        # TODO: handle situations with list-valued id types
    df.drop_duplicates(inplace=True)
    df = df.groupby(id_column).mean()
    # df.set_index(id_column, inplace=True)


    scores = {ID: float(df.ix[ID, score_column]) for ID in df.index}
    nan_keys = {ID for ID in scores if math.isnan(scores[ID])}
    for ID in nan_keys:
        del scores[ID]
    for ID in scores.keys():
        # handle entrez IDs and other integer like IDs
        try:
            _ID = str(int(float(ID)))
            if _ID != ID:
                scores[_ID] = scores[ID]
                del scores[ID]
        except:
            pass
    return scores

def parse_geneset(layer_name, argparse_args, id_mapper):
    # please do not use this to do funny things
    return _parse_geneset(eval('argparse_args.'+layer_name+'_file'),
                          eval('argparse_args.'+layer_name+'_genesets'),
                          eval('argparse_args.'+layer_name),
                          eval('argparse_args.'+layer_name+'_id_type'),
                          argparse_args.graph_id_type,
                          id_mapper)

def _parse_geneset(geneset_files, geneset_names, genes, geneset_id_type, target_id_type, id_mapper):
    if geneset_files is None and genes is None:
        return set()
    S = geneset.GeneSet('S', set(genes.split(',')))
    geneset_names_cnt = 0
    for geneset_file in geneset_files:
        suffix = geneset_file.split('.')[-1]
        if suffix != 'gmt':
            S |= geneset.GeneSet.from_grp(geneset_file)
        else:
            coll = geneset.GeneSetCollection.from_gmt(geneset_file)
            names = geneset_names[geneset_names_cnt].split(',')
            S |= coll.union(names)
            geneset_names_cnt += 1
    if geneset_id_type != graph_id_type:
        S.map(id_mapper, geneset_id_type, graph_id_type, inplace=True)
    return S
