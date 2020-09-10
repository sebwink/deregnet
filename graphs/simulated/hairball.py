import sys
import random
import networkx as nx

COLORS = ['RED', 'BLUE']

def main(out):
    G = nx.erdos_renyi_graph(200, 0.2)
    for node in G.nodes:
        G.nodes[node]['color'] = random.choice(COLORS)
    nx.write_graphml(G, out)

if __name__ == '__main__':
    main(sys.argv[1])
