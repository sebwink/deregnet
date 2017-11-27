import numpy as np
import igraph as ig

def shortest_path_distance_score(method):

    def _method(self, subgraph1, subgraph2, match_on, score_attr, **kwargs):
        mode = kwargs.get('mode', ig.ALL)
        if self.d[mode] is None:
            self.d[mode] = self.calculate_shortest_path_distances(**kwargs)
        return method(subgraph1, subgraph2, match_on, score_attr, **kwargs)

    return _method


class SubgraphProximityScorer:
    '''
    Class implementing some concepts for quantifying "proximity"
    of two subgraphs in the context of their base graph.

    Most efficient for batch scoring, i.e. scoring as many subgraph
    pairs with one instance of this class. (Every time it computes
    a score based on shortest path distances for the first time
    shortest path distances for every vertex pair in the base graph
    is computed.) "Most efficient" can also mean it is very in-
    efficient, for example in the case when all scored subgraph
    pairs are 'close by' and one could have come by with much less
    shortest path calculations than calculate them globally for all
    vertex pairs in the base graph.

    Also, all implementations below involving a vertex score weighting
    assume that the high vertex scores mean higher deregulation. While
    often true, e.g. fold-change, differential expression indicator, etc.
    scores with opposite semantics are also possible, i.e. scores based
    on p-values.
    '''
    # TODO: implement adaptive shortest path distance calculation based
    # the actual vertex sets of the scored subgraphs and cache these 
    # distances; the same for whole paths
    # TODO: support score weighting for vertex scores with minimization 
    # semantics

    def __init__(self, base_graph):
        self.graph = base_graph
        # 1 = ig.OUT, 2 = ig.IN, 3 = ig.ALL
        self.d = { 1: None, 2: None, 3: None }

    @classmethod
    def vertex_intersection_score(cls, subgraph1, subgraph2, match_on='name'):
        '''     _
        |V(G1) | | V(G2)|

        Args:

            subgraph1 (ig.Graph): first subgraph
            subgraph2 (ig.Graph): second subgraph
            match_on (str): vertex attribute to match vertices on

        Return:

            int: vertex intersection score
        '''
        return len(  { v[match_on] for v in subgraph1.vs } \
                   & { v[match_on] for v in subgraph2.vs } )

    @shortest_path_distance_score
    def unweighted_shortest_path_score(self,
                                       subgraph1,
                                       subgraph2,
                                       match_on='name',
                                       score_attr='score',
                                       **kwargs):
        '''
                                    _                _
                     ---           |                  |
         1     1     \             |    1        1    |
        --- -------- /   f(u,v) *  | ------- + ------ |
         2  N(G1,G2) ---           |  d(u,v)   d(v,u) |
                     u,v           |_                _|
                     u!=v

        with N(G1,G2) = |V(G1)| + |V(G2)| - |V(G1) | | V(G2)|

        and f(u,v) = 1

        In case of mode = ig.ALL we have d(u,v) = d(v,u) this reduces to:

                 ---
           1     \     f(u,v)
        -------- /    -------
        N(G1,G2) ---   d(u,v)
                 u,v
                 u!=v

        Args:

            subgraph1 (ig.Graph): first subgraph
            subgraph2 (ig.Graph): second subgraph
            match_on (str): vertex attribute to match vertices on

        Return:

            int: unweighted shorted path score
        '''
        pass

    def uwsps(self, subgraph1, subgraph2, **kwargs):
        pass

    @shortest_path_distance_score
    def endpoint_weighted_shortest_path_score(self,
                                              subgraph1,
                                              subgraph2,
                                              match_on='name',
                                              score_attr='score',
                                              **kwargs):
        '''
                                    _                _
                     ---           |                  |
         1     1     \             |    1        1    |
        --- -------- /   f(u,v) *  | ------- + ------ |
         2  N(G1,G2) ---           |  d(u,v)   d(v,u) |
                     u,v           |_                _|
                     u!=v

        with N(G1,G2) = |V(G1)| + |V(G2)| - |V(G1) | | V(G2)|

                      s(u) + s(v)
        and f(u,v) = -------------
                           2

        In case of mode = ig.ALL we have d(u,v) = d(v,u) this reduces to:

                 ---
           1     \     f(u,v)
        -------- /    -------
        N(G1,G2) ---   d(u,v)
                 u,v
                 u!=v

        Args:

            subgraph1 (ig.Graph): first subgraph
            subgraph2 (ig.Graph): second subgraph
            match_on (str): vertex attribute to match vertices on

        Return:

            int: endpoint weighted shorted path score

        '''
        pass

    def ewsps(self, subgraph1, subgraph2, **kwargs):
        '''
        alias for endpoint_weighted_shortest_path_score
        '''
        return self.endpoint_weighted_shortest_path_score(subgraph1, subgraph2, **kwargs)

    @shortest_path_score
    def path_weighted_shortest_path_score(self,
                                          subgraph1,
                                          subgraph2,
                                          match_on='name',
                                          score_attr='score',
                                          **kwargs):
        '''

                     ---
         1     1     \                1                  1
        --- -------- /   f1(u,v) * ------- + f2(u,v) * ------
         2  N(G1,G2) ---            d(u,v)             d(v,u)
                     u,v
                     u!=v

        with N(G1,G2) = |V(G1)| + |V(G2)| - |V(G1) | | V(G2)|

                       ---                      ---
                       \                        \
        and f1(u,v) =  /    s(w)  , f2(u,v) =   /     s(w)
                       ---                      ---
                   w in p(u,v)              w in p(v,u)

        where p(x,y) is the shortest path x --> y.

        In case of mode = ig.ALL we have d(u,v) = d(v,u) this reduces to:

                 ---
           1     \     f(u,v)
        -------- /    -------
        N(G1,G2) ---   d(u,v)
                 u,v
                 u!=v

        Args:

            subgraph1 (ig.Graph): first subgraph
            subgraph2 (ig.Graph): second subgraph
            match_on (str): vertex attribute to match vertices on

        Return:

            int: path weighted shorted path score
        '''
        pass

    def pwsps(self, subgraph1, subgraph2, **kwargs):
        pass

    def calculate_shortest_path_distances(self, **kwargs):
        return np.array(base_graph.shortest_paths(self.graph.vs, self.graph.vs, **kwargs))

    def _init_d_all(self, **kwargs):
        self.d_all = self.calculate_shortest_path_distances(mode=ig.ALL, **kwargs)

    def _init_d_out(self, **kwargs):
        self.d_out = self.calculate_shortest_path_distances(mode=ig.OUT, **kwargs)

    def _init_d_in(self, **kwargs):
        self.d_in = self.calculate_shortest_path_distances(mode=ig.IN, **kwargs)

    @classmethod
    def vertex_union_size(cls, subgraph1, subgraph2, match_on='name'):
        return len(  { v[match_on] for v in subgraph1.vs } \
                   | { v[match_on] for v in subgraph2.vs } )
