import requests
import igraph as ig

OMNIPATH_ROOT='http://omnipathdb.org'

def get_omnipath_interaction_graph():
    query = OMNIPATH_ROOT+'/interactions'
    response = requests.get(query)
    if response.status_code != 200:
        print('Error requesting %s ...' % query)
        return None
    data = response.content.decode('utf-8')
    data = [line.strip() for line in data.split('\n') if line]
    nodes = set()
    edges = list()
    is_stimulation = list()
    is_inhibition = list()
    for line in data[1:]:
        source, target, is_directed, _is_stimulation, _is_inhibition = line.split('\t')[:5]
        is_directed = bool(int(is_directed))
        if not is_directed:
            continue
        print('%s -->  %s' % (source, target))
        nodes.add(source)
        nodes.add(target)
        _is_stimulation = bool(int(_is_stimulation))
        _is_inhibition = bool(int(_is_inhibition))
        edges.append((source, target))
        is_stimulation.append(_is_stimulation)
        is_inhibition.append(_is_inhibition)
    nodes = list(nodes)
    nodes2index = { node : nodes.index(node) for node in nodes }
    edges = [(nodes2index[edge[0]], nodes2index[edge[1]]) for edge in edges]
    graph = ig.Graph(directed=True)
    graph.add_vertices(len(nodes))
    graph.vs['name'] = nodes
    graph.add_edges(edges)
    graph.es['is_stimulation'] = is_stimulation
    graph.es['is_inhibition'] = is_inhibition
    return graph

if __name__ == '__main__':
    graph = get_omnipath_interaction_graph()
    graph.write_graphml('../omnipath_directed_interactions.graphml')
