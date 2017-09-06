import requests
import igraph as ig

OMNIPATH_ROOT='http://omnipathdb.org'

def get_omnipath_ptm_graph():
    query = OMNIPATH_ROOT+'/ptms'
    response = requests.get(query)
    if response.status_code != 200:
        print('Error requesting %s ...' % query)
        return None
    data = response.content.decode('utf-8')
    data = [line.strip() for line in data.split('\n') if line]
    nodes = set()
    edges = list()
    residue_type = list()
    residue_offset = list()
    modification = list()
    for line in data[1:]:
        source, target, _residue_type, _residue_offset, _modification = line.split('\t')[:5]
        print('%s --> %s, %s' % (source, target, _modification))
        nodes.add(source)
        nodes.add(target)
        edges.append((source, target))
        residue_type.append(_residue_type)
        residue_offset.append(_residue_offset)
        modification.append(_modification)
    nodes = list(nodes)
    nodes2index = { node : nodes.index(node) for node in nodes }
    edges = [(nodes2index[edge[0]], nodes2index[edge[1]]) for edge in edges]
    graph = ig.Graph(directed=True)
    graph.add_vertices(len(nodes))
    graph.vs['name'] = nodes
    graph.add_edges(edges)
    graph.es['residue_type'] = residue_type
    graph.es['residue_offset'] = residue_offset
    graph.es['modification'] = modification
    return graph

if __name__ == '__main__':
    graph = get_omnipath_ptm_graph()
    graph.write_graphml('../omnipath_ptm_graph.graphml')
