#include <set>
#include <utility>
#include <functional>

#include <deregnet/usinglemon.h>
#include <deregnet/DeregnetStartHeuristic.h>

using namespace std;

namespace deregnet {

DeregnetStartHeuristic::DeregnetStartHeuristic(Graph* xgraph,
                                               NodeMap<double>* xscore,
                                               Node* xroot,
                                               int xsize,
                                               set<Node>* xexclude,
                                               set<Node>* xreceptors,
                                               function<bool(double, double)> xcmp)
 : graph { xgraph },
   score { xscore },
   root { xroot },
   size { xsize },
   exclude { xexclude },
   receptors { xreceptors },
   cmp { xcmp },
   start_solution { new set<Node> }
{  }

bool DeregnetStartHeuristic::run() {
    if (!root)
        root = get_best_root();
    if (!root)
        return false;
    found_node = true;
    start_solution->insert(*root);
    Node* next { nullptr };
    while ( search_further() ) {
        next = get_next_node();
        if (  next ) {
            start_solution->insert(*next);
            found_node = true;
        }
        else
            found_node = false;
    }
    return success;
}

pair<Node, set<Node>>* DeregnetStartHeuristic::getStartSolution() {
    return new pair<Node, set<Node>>( make_pair(*root, *start_solution) );
}

Node* DeregnetStartHeuristic::get_best_root() {
    double best;
    Node* best_node { nullptr };
    if (receptors)
        for (auto v : *receptors)
            update_best_node(&best_node, &v, &best);
    else
        for (NodeIt v(*graph); v != INVALID; ++v)
            update_best_node(&best_node, &v, &best);
    return best_node;
}

Node* DeregnetStartHeuristic::get_next_node() {
    double best;
    Node* best_node { nullptr };
    Node* u = new Node();
    for (auto v : *start_solution) {
        for (OutArcIt a(*graph, v); a != INVALID; ++a) {
            *u = graph->target(a);
            if ( feasible_node(u) )
                update_best_node(&best_node, u, &best);
        }
    }
    return best_node;
}

void DeregnetStartHeuristic::update_best_node(Node** current_best_node, Node* node, double* best) {
    if (!(*current_best_node)) {
        *current_best_node = new Node();
        **current_best_node = *node;
        *best = (*score)[*node];
    }
    else if ( cmp(*best, (*score)[*node]) ) {
        **current_best_node = *node;
        *best = (*score)[*node];
    }
}

}    //    namespace deregnet
