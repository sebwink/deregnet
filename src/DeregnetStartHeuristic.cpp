// --------------------------------------------------------------------------
//                   deregnet -- find deregulated pathways
// --------------------------------------------------------------------------
// Copyright Sebastian Winkler --- Eberhard Karls University Tuebingen, 2016
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Sebastian Winkler $
// $Authors: Sebastian Winkler $
// --------------------------------------------------------------------------
//

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
                                               set<Node>* xexclude,
                                               set<Node>* xreceptors,
                                               function<bool(double, double)> xcmp)
 : graph { xgraph },
   score { xscore },
   root { nullptr },
   original_root { xroot },
   exclude { xexclude },
   receptors { xreceptors },
   cmp { xcmp },
   start_solution { new set<Node> }
{
    if ( xroot ) {
        root = new Node();
        *root = *xroot;
    }
}

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
        if ( next ) {
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
