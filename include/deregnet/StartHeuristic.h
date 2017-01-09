#ifndef DEREGNET_START_HEURISTIC_H
#define DEREGNET_START_HEURISTIC_H

#include <set>

#include <deregnet/usinglemon.h>

class StartHeuristic {

  private:

    Graph* graph;
    NodeMap<double>* score;
    int iterations;
    std::set<Node> initial_solution;

  public:

    StartHeuristic(Graph* xgraph, NodeMap<double>* xscore, int iterations);
    void run();
    std::set<Node> getInitialSolution();


  private:

    void getNextMaximalNode(Node* max_node, double max_score);

};

#endif    //    DEREGNET_START_HEURISTIC
