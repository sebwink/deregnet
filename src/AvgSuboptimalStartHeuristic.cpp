// --------------------------------------------------------------------------
//                   deregnet -- find deregulated pathways
// --------------------------------------------------------------------------
// Copyright Sebastian Winkler --- Eberhard Karls University Tuebingen, 2020
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

#include "usinglemon.hpp"
#include "AvgSuboptimalStartHeuristic.hpp"

using namespace std;

namespace deregnet {

AvgSuboptimalStartHeuristic::AvgSuboptimalStartHeuristic(Graph* xgraph,
                                                         NodeMap<double>* xscore,
                                                         Node* root,
                                                         std::set<Node>* exclude,
                                                         std::set<Node>* receptors,
                                                         std::function<bool(double, double)> xcmp,
                                                         NodeMap<std::string>* xnodeid,
                                                         std::set<std::string>* xnodes_sof_far,
                                                         double xmax_overlap,
                                                         int xmin_size,
                                                         int xmax_size)
 : DeregnetStartHeuristic(xgraph, xscore, root, exclude, receptors, xcmp),
   nodeid { xnodeid },
   nodes_so_far { xnodes_sof_far },
   max_overlap { xmax_overlap },
   min_size { xmin_size },
   max_size { xmax_size }
{ }

Node* AvgSuboptimalStartHeuristic::get_best_root() {
    double best;
    Node* best_node { nullptr };
    bool overlap_matters { max_overlap < 1.0 / nodes_so_far->size() };
    if (receptors)
        for (auto v : *receptors) {
            if ( is_overlap_node(&v) && overlap_matters )
                continue;
            else
                update_best_node(&best_node, &v, &best);
        }
    else
        for (NodeIt v(*graph); v != INVALID; ++v) {
            if ( is_overlap_node(&v) && overlap_matters )
                continue;
            else
                update_best_node(&best_node, &v, &best);
        }
    return best_node;
}

bool AvgSuboptimalStartHeuristic::search_further() {
    int current_size { static_cast<int>(start_solution->size()) };
    if (current_size < min_size) {
        if (found_node)
            return true;
        else
            return false;
    }
    if (current_size >= min_size && current_size < max_size) {
        if (found_node)
            return true;
        else {
            success = true;
            return false;
        }
    }
    if (current_size == max_size) {
        success = true;
        return false;
    }
    return false;
}

bool AvgSuboptimalStartHeuristic::feasible_node(Node* node) {
    double new_avg { update_avg(node) };
    int current_size { static_cast<int>(start_solution->size()) };
    if (new_avg < current_avg && current_size > min_size)
        return false;
    current_avg = new_avg;

    bool overlap_node { is_overlap_node(node) };
    double potential_overlap { static_cast<double>( current_overlapping_nodes / start_solution->size() ) };
    if (exclude) {
        if (overlap_node) {
            if (   start_solution->find(*node) == start_solution->end()
                && potential_overlap < max_overlap
                && exclude->find(*node) == exclude->end()) {
                current_overlapping_nodes++;
                return true;
            }
        }
        else {
            if (start_solution->find(*node) == start_solution->end() && exclude->find(*node) == exclude->end())
                return true;
        }
    }
    else {
        if (overlap_node) {
            if (   start_solution->find(*node) == start_solution->end()
                && potential_overlap < max_overlap) {
                current_overlapping_nodes++;
                return true;
            }
        }
        else {
            if (start_solution->find(*node) == start_solution->end())
                return true;
        }
    }
    return false;
}

bool AvgSuboptimalStartHeuristic::is_overlap_node(Node* node) {
    for (set<string>::iterator nodeID = nodes_so_far->begin(); nodeID != nodes_so_far->end(); ++nodeID)
        if ( *nodeID == (*nodeid)[*node] )
            return true;
    return false;
}

double AvgSuboptimalStartHeuristic::update_avg(Node* node) {
    int current_size { static_cast<int>(start_solution->size()) };
    return ((*score)[*node] + current_size * current_avg) / (current_size + 1);
}

}    //    namespace deregnet
