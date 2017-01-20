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

#include <cmath>
#include <iostream>
#include <map>
#include <set>
#include <utility>

#include <gurobi_c++.h>

#include <lemon/connectivity.h>
#include <lemon/adaptors.h>

#include <deregnet/usinglemon.h>

#include <deregnet/LazyConstraintCallback.h>

using namespace std;

namespace deregnet {

LazyConstraintCallback::LazyConstraintCallback(std::map<Node, GRBVar>* xx,
                                               Graph* xgraph,
                                               Node* xroot,
                                               double* xgap_cut)
 : x { xx },
   graph { xgraph },
   original_root { xroot },
   root { nullptr },
   gap_cut { xgap_cut }
 {
    if (xroot) {
        root = new Node();
        *root = *original_root;
    }
 }

void LazyConstraintCallback::callback() {
    try {
        if (where == GRB_CB_MIPSOL) { // new incumbent found
            NodeFilter solution_filter(*graph);
            selected_nodes = {};
            get_solution_filter(solution_filter);
            InducedSubgraph current_subgraph(*graph, solution_filter);
            InducedSubgraph::NodeMap<int> component_map(current_subgraph);
            int num_components { lemon::stronglyConnectedComponents(current_subgraph, component_map) };
            if (num_components > 1)
                check_and_set_lazy_constr(num_components, component_map);
        }
        else if (where == GRB_CB_MIP && gap_cut) {
            double best_objective { getDoubleInfo(GRB_CB_MIP_OBJBST) };
            double best_bound { getDoubleInfo(GRB_CB_MIP_OBJBND) };
            if ( abs(best_objective - best_bound) < (*gap_cut) * (1.0 + abs(best_objective)) ) {
                cout << "Achieved " << *gap_cut << " gap. Stopping optimization." << endl;
                abort();
            }
        }
    }
    catch (GRBException e) {
        std::cout << "Error number: " << e.getErrorCode() << std::endl;
        std::cout << e.getMessage() << std::endl;
    }
    catch (...) {
        std::cout << "Error during callback" << std::endl;
    }
}

void LazyConstraintCallback::check_and_set_lazy_constr(const int num_components,
                                                       const InducedSubgraph::NodeMap<int>& component_map) {
    for(int k = 0; k < num_components; k++) {
        std::set<Node> component;
        std::set<Node> parents;
        std::set<Node> global_parents;
        bool more_than_one_node = get_component_nodes(component_map, component, k);
        if (!is_root_component(component) && more_than_one_node) {
            get_parents(component, parents, global_parents);
            if (parents.size() == 0)
                set_lazy_constraint(component, global_parents);
        }
    }
}

bool LazyConstraintCallback::get_component_nodes(const InducedSubgraph::NodeMap<int>& component_map,
                                                 std::set<Node>& component,
                                                 const int k) {
    for (auto v : selected_nodes)
        if (k == component_map[v])
            component.insert(v);
    if (component.size() > 1)
        return true;
    return false;
}

bool LazyConstraintCallback::is_root_component(const std::set<Node>& component) {
    for (auto v : component)
        if (v == *root)
            return true;
    return false;
}

void LazyConstraintCallback::get_parents(const std::set<Node>& component,
                                         std::set<Node>& parents,
                                         std::set<Node>& global_parents) {
    for (auto v : component) {
        for (InArcIt a(*graph, v); a != INVALID; ++a) {
            Node u = graph->source(a);
            if (component.find(u) == component.end()) {
                global_parents.insert(u);
                if (find(selected_nodes.begin(), selected_nodes.end(), u) != selected_nodes.end())
                    parents.insert(u);
            }
        }
    }
}


// LazyConstraintCallbackRoot implementation

void LazyConstraintCallbackRoot::get_solution_filter(NodeFilter& solution_filter) {
    for (NodeIt v(*graph); v != INVALID; ++v) {
        if (getSolution((*x)[v]) >= 0.98) {
            solution_filter[v] = true;
            selected_nodes.push_back(v);
        }
        else
            solution_filter[v] = false;
    }
}

void LazyConstraintCallbackRoot::set_lazy_constraint(const std::set<Node>& component,
                                                     const std::set<Node>& global_parents) {
    GRBLinExpr component_expr;
    GRBLinExpr global_parents_expr;
    for (auto v : component)
        component_expr += (*x)[v];
    for (auto v : global_parents)
        global_parents_expr += (*x)[v];
    GRBLinExpr lazy_constr_lhs = component_expr - global_parents_expr;
    int comp_size = component.size();
    addLazy(lazy_constr_lhs <= comp_size - 1);
}


// LazyConstraintCallbackNoRoot implementation

LazyConstraintCallbackNoRoot::LazyConstraintCallbackNoRoot(std::map<Node, GRBVar>* xx,
                                                           std::map<Node, GRBVar>* yy,
                                                           Graph* xgraph,
                                                           double* gap_cut)
 : LazyConstraintCallback(xx, xgraph, nullptr, gap_cut),
   y { yy }
{ }



void LazyConstraintCallbackNoRoot::get_solution_filter(NodeFilter& solution_filter) {
    for (NodeIt v(*graph); v != INVALID; ++v) {
        if (getSolution((*x)[v]) >= 0.98) {
            solution_filter[v] = true;
            selected_nodes.push_back(v);
            if (getSolution((*y)[v]) >= 0.98)
                root = &selected_nodes.back();
        }
        else
            solution_filter[v] = false;
    }
}


void LazyConstraintCallbackNoRoot::set_lazy_constraint(const std::set<Node>& component,
                                                       const std::set<Node>& global_parents) {
    GRBLinExpr component_expr;
    GRBLinExpr global_parents_expr;
    for (auto v : component)                      // void build_component_expr(GRBLinExpr& expr); ...
        component_expr += (*x)[v] - (*y)[v];
    for (auto v : global_parents)
        global_parents_expr += (*x)[v];
    GRBLinExpr lazy_constr_lhs = component_expr - global_parents_expr;
    int comp_size = component.size();
    addLazy(lazy_constr_lhs <= comp_size - 1);
}

}   // namespace deregnet
