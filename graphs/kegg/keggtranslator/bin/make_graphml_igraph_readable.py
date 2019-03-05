import sys
import networkx as nx

def main(graphml):
    g = nx.read_graphml(graphml)
    nx.write_graphml(g, graphml)

if __name__ == '__main__':
    main(sys.argv[1])
