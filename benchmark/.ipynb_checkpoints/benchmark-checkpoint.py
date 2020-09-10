import os
import time
import random
import json
import igraph as ig

from deregnet.core import SubgraphFinder

SELF = os.path.dirname(__file__)
GRAPHML = os.getenv('GRAPHML', os.path.join(SELF, 'data/kegg_hsa.graphml'))

MIN_SUBGRAPH_SIZE = os.getenv('MIN_SUBGRAPH_SIZE', 5)
MAX_SUBGRAPH_SIZE = os.getenv('MAX_SUBGRAPH_SIZE', 15)

ACTUAL_MIN_SUBGRAPH_SIZE = os.getenv(
    'ACTUAL_MIN_SUBGRAPH_SIZE',
    5,
)
ACTUAL_MAX_SUBGRAPH_SIZE = os.getenv(
    'ACTUAL_MAX_SUBGRAPH_SIZE',
    15,
)

FLIP_P = os.getenv('FLIP_P', 0.01)

NRUNS = os.getenv('RUNS', 100)
TIME_LIMIT_PER_RUN = os.getenv('TIME_LIMIT_PER_RUN', 600)

RESULTS = os.getenv(
    'RESULTS',
    'benchmark.%d.%d.%d.%d.%s.%d.%d.json' % (
        int(MIN_SUBGRAPH_SIZE),
        int(ACTUAL_MIN_SUBGRAPH_SIZE),
        int(MAX_SUBGRAPH_SIZE),
        int(ACTUAL_MAX_SUBGRAPH_SIZE),
        str(FLIP_P),
        int(NRUNS),
        int(TIME_LIMIT_PER_RUN),
    ),
)


LOGS = os.getenv(
    'LOGS',
    'benchmark.%d.%d.%d.%d.%s.%d.%d.log' % (
        int(MIN_SUBGRAPH_SIZE),
        int(ACTUAL_MIN_SUBGRAPH_SIZE),
        int(MAX_SUBGRAPH_SIZE),
        int(ACTUAL_MAX_SUBGRAPH_SIZE),
        str(FLIP_P),
        int(NRUNS),
        int(TIME_LIMIT_PER_RUN),
    ),
)


def simulate_subgraph(graph):
    num_nodes = len(graph.vs)
    root_index = random.randint(0, num_nodes)
    subgraph_size = random.randint(
        ACTUAL_MIN_SUBGRAPH_SIZE,
        ACTUAL_MAX_SUBGRAPH_SIZE,
    )

    MAX_ITER = 4 * ACTUAL_MAX_SUBGRAPH_SIZE

    def add_node(subgraph):
        if len(subgraph) > subgraph_size:
            return False
        randind = random.randrange(len(subgraph))
        node = list(graph.vs.select(name_eq=subgraph[randind]))
        if not node:
            return True
        node = node[0]
        for n in (graph.vs['name'][n] for n in graph.neighbors(node, ig.OUT)):
            if n not in subgraph:
                subgraph.append(n)
                return True

    subgraph = [graph.vs[root_index]['name']]
    for iters in range(MAX_ITER):
        if not add_node(subgraph):
            break
    if len(subgraph) < ACTUAL_MIN_SUBGRAPH_SIZE:
        return None
    return subgraph


def simulate_scores(subgraph, graph):
    def flip(p, score):
        u = random.random()
        return score if u > p else abs(score - 1)
    scores = {v: 1 if v in subgraph else 0 for v in graph.vs['name']}
    return {v: flip(FLIP_P, score) for v, score in scores.items()}


def run_drgnt(finder, scores):
    def union_subgraph(subgraphs):
        if not subgraphs:
            return None
        seed = set(subgraphs[0].vs['name'])
        if len(subgraphs) > 1:
            return seed.union(
                *[
                    set(sg.vs['name'])
                    for sg in subgraphs[2:]
                ]
            )
        else:
            return seed
    start = time.time()
    remaining_time = TIME_LIMIT_PER_RUN
    subgraphs = []
    for size in range(MIN_SUBGRAPH_SIZE, MAX_SUBGRAPH_SIZE+1):
        try:
            result = finder.run_absolute_deregnet(
                scores=scores,
                size=size,
                time_limit=remaining_time,
            )
            subgraphs.append(result.subgraphs[0])
            remaining_time = TIME_LIMIT_PER_RUN - (time.time() - start)
        except Exception:
            continue
    return union_subgraph(subgraphs)


def run_avgdrgnt_gcc(finder, scores):
    try:
        result = finder.run_average_deregnet(
            scores=scores,
            min_size=MIN_SUBGRAPH_SIZE,
            max_size=MAX_SUBGRAPH_SIZE,
            algorithm='gcc',
            time_limit=TIME_LIMIT_PER_RUN,
        )
        return set(result.subgraphs[0].vs['name'])
    except Exception:
        return None


def run_avgdrgnt_dta(finder, scores):
    try:
        result = finder.run_average_deregnet(
            scores=scores,
            min_size=MIN_SUBGRAPH_SIZE,
            max_size=MAX_SUBGRAPH_SIZE,
            algorithm='dta',
            time_limit=TIME_LIMIT_PER_RUN / 3,  # TODO: not optimal
        )
        return set(result.subgraphs[0].vs['name'])
    except Exception:
        return None


def benchmark(runner, subgraph, finder, scores):
    data = {}
    start = time.time()
    result = runner(finder, scores)
    data['subgraph_size'] = len(subgraph)
    if not result:
        data['success'] = False
        data['hits'] = None
        data['fp'] = None
    else:
        data['success'] = True
        data['hits'] = len(subgraph.intersection(result))
        data['fp'] = len(result.difference(subgraph))
    data['time'] = time.time() - start
    return data


if __name__ == '__main__':
    graph = ig.Graph.Read_GraphML(GRAPHML)
    data = {
        'known_min_size': MIN_SUBGRAPH_SIZE,
        'min_size': ACTUAL_MIN_SUBGRAPH_SIZE,
        'known_max_size': MAX_SUBGRAPH_SIZE,
        'max_size': ACTUAL_MAX_SUBGRAPH_SIZE,
        'flip_p': FLIP_P,
        'nruns': NRUNS,
        'time_limit': TIME_LIMIT_PER_RUN,
        'avg_gcc': [],
        'avg_dta': [],
        'abs': [],
        'subgraph': [],
        'score': [],
    }
    nruns = 0
    while nruns < NRUNS:
        subgraph = simulate_subgraph(graph)
        if subgraph is None:
            continue
        nruns += 1
        data['subgraph'].append(subgraph)
        subgraph = set(subgraph)
        scores = simulate_scores(subgraph, graph)
        data['score'].append(scores)
        finder = SubgraphFinder(graph, log_file=LOGS)
        data['abs'].append(
            benchmark(run_drgnt, subgraph, finder, scores),
        )
        data['avg_gcc'].append(
            benchmark(run_avgdrgnt_gcc, subgraph, finder, scores),
        )
        data['avg_dta'].append(
            benchmark(run_avgdrgnt_dta, subgraph, finder, scores),
        )
    with open(RESULTS, 'w') as fp:
        json.dump(data, fp)
