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

#ifndef DEREGNET_MODEL_H
#define DEREGNET_MODEL_H

#include <iostream>
#include <string>
#include <map>
#include <set>
#include <utility>
#include <vector>

#include <lemon/connectivity.h>
#include <lemon/adaptors.h>

#include <gurobi_c++.h>

#include <grbfrc/FMILP.h>

#include <deregnet/utils.h>
#include <deregnet/usinglemon.h>
#include <deregnet/LazyConstraintCallback.h>

namespace deregnet {

template <typename Model>
class DeregnetModel {

  protected:

    GRBEnv env;
    Model model { Model(env) };
    std::map<Node, GRBVar> x;
    Graph* graph;
    NodeMap<double>* score;
    NodeMap<std::string>* nodeid;
    Node* root { nullptr };
    std::map<Node, GRBVar>* y { nullptr };

  public:

    DeregnetModel(Graph* graphp,
                  NodeMap<double>* scorep,
                  NodeMap<string>* nodeidp,
                  Node* rootp);
    void createVariables();
    void addBaseConstraints(int size);
    void addIncludeConstraints(std::set<Node>* include);
    void addExcludeConstraints(std::set<Node>* exclude);
    void addReceptorConstraints(std::set<Node>* receptors);
    void addTerminalConstraints(std::set<Node>* terminals);
    bool solve(std::pair<Node, std::set<Node>>* start_solution,
               double* time_limit,
               double* gap_cut,
               std::string model_sense,
               grbfrc::Algorithm algorithm = grbfrc::Algorithm::GCC);
    void addSuboptimalityConstraint(std::set<std::string> nodes_so_far,
                                    double max_overlap,
                                    int size);
    Solution getCurrentSolution();

  private:

    void createVariablesRoot();                 // specialize for GRBModel and FMILP
    void createVariablesNoRoot();               // specialize for GRBModel and FMILP
    void addBaseConstraintsRoot(int size);
    void addBaseConstraintsNoRoot(int size);
    void setStartSolution(std::pair<Node, std::set<Node>>* start_solution);  // specialize for GRBModel and FMILP ?
    void solve_common(std::pair<Node, std::set<Node>>* start_solution,
                      double* time_limit,
                      double* gap_cut,
                      std::string model_sense);
    void setCallbackRoot(double* gap_cut);                     // specialize for GRBModel and FMILP ?!
    void setCallbackNoRoot(double* gap_cut);                   // specialize for GRBModel and FMILP ?!

};

template <typename Model>
DeregnetModel<Model>::DeregnetModel(Graph* graphp,
                                    NodeMap<double>* scorep,
                                    NodeMap<std::string>* nodeidp,
                                    Node* rootp)
 : graph { graphp },
   score { scorep },
   nodeid { nodeidp },
   root { rootp }
{ }

template <typename Model>
void DeregnetModel<Model>::createVariables() {
    if (root)
        createVariablesRoot();
    else
        createVariablesNoRoot();
    model.update();
}

template <typename Model>
void DeregnetModel<Model>::addBaseConstraints(int size) {
    if (root)
        addBaseConstraintsRoot(size);
    else
        addBaseConstraintsNoRoot(size);
    model.update();
}

template <typename Model>
void DeregnetModel<Model>::addBaseConstraintsRoot(int size) {
    GRBLinExpr subgraph_size_lhs;
    for (NodeIt v(*graph); v != INVALID; ++v) {
        subgraph_size_lhs += x[v];
        GRBLinExpr parent_sum;
        for (InArcIt a(*graph,v); a != INVALID; ++a)
            parent_sum += x[graph->source(a)];
        if (v != *root)
            model.addConstr(x[v] - parent_sum <= 0);
    }
    model.addConstr(subgraph_size_lhs == size);
    model.addConstr(x[*root] == 1);
}

template <typename Model>
void DeregnetModel<Model>::addBaseConstraintsNoRoot(int size) {
    GRBLinExpr subgraph_size_lhs;
    GRBLinExpr exactly_one_root_lhs;
    for (NodeIt v(*graph); v != INVALID; ++v) {
        model.addConstr((*y)[v] <= x[v]);
        subgraph_size_lhs += x[v];
        exactly_one_root_lhs += (*y)[v];
        GRBLinExpr parent_sum;
        for (InArcIt a(*graph,v); a != INVALID; ++a)
            parent_sum += x[graph->source(a)];
        model.addConstr(x[v] - (*y)[v] - parent_sum <= 0);
    }
    model.addConstr(subgraph_size_lhs == size);
    model.addConstr(exactly_one_root_lhs == 1);
}

template <typename Model>    // remains to be seen if we need a specialization for FMILP
void DeregnetModel<Model>::setStartSolution(std::pair<Node, std::set<Node>>* start_solution) {
    for (NodeIt v(*graph); v != INVALID; ++v) {
        if ((start_solution->second).find(v) != (start_solution->second).end())
            x[v].set(GRB_DoubleAttr_Start, 1.0);
        else
            x[v].set(GRB_DoubleAttr_Start, 0.0);
        if (!root) {
            if (v == start_solution->first)
                (*y)[v].set(GRB_DoubleAttr_Start, 1.0);
            else
                (*y)[v].set(GRB_DoubleAttr_Start, 0.0);
        }
    }
}

template <typename Model>
void DeregnetModel<Model>::addIncludeConstraints(std::set<Node>* include) {
    for (Node v : *include)
        model.addConstr(x[v] == 1);
    model.update();
}

template <typename Model>
void DeregnetModel<Model>::addExcludeConstraints(std::set<Node>* exclude) {
    for (Node v : *exclude)
        model.addConstr(x[v] == 0);
    model.update();
}

template <typename Model>
void DeregnetModel<Model>::addReceptorConstraints(std::set<Node>* receptors) {
    if (!root)
        for (NodeIt v(*graph); v != INVALID; ++v)
            if (receptors->find(v) == receptors->end())
                model.addConstr((*y)[v] == 0);
    model.update();
}

template <typename Model>
void DeregnetModel<Model>::addTerminalConstraints(std::set<Node>* terminals) { // consider using complement of terminals ...
    for (NodeIt v(*graph); v != INVALID; ++v)
        if (terminals->find(v) == terminals->end()) {
            GRBLinExpr out_neighbor_expr;
            for (OutArcIt a(*graph, v); a != INVALID; ++a)
                out_neighbor_expr += x[graph->target(a)];
            model.addConstr(x[v] - out_neighbor_expr <= 0);
        }
}

template <typename Model>
void DeregnetModel<Model>::solve_common(std::pair<Node, std::set<Node>>* start_solution,
                                        double* time_limit,
                                        double* gap_cut,
                                        std::string model_sense) {
    if (model_sense == "min")
        model.set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);
    else
        model.set(GRB_IntAttr_ModelSense, GRB_MAXIMIZE);

    model.set(GRB_IntParam_LazyConstraints, 1);

    if (time_limit)
        model.set(GRB_DoubleParam_TimeLimit, *time_limit);

    if (root) {
        if (start_solution)
            setStartSolution(start_solution);
        setCallbackRoot(gap_cut);
    }
    else {
        if (start_solution)
            setStartSolution(start_solution);
        setCallbackNoRoot(gap_cut);
    }
    model.update();
}

template <typename Model>
void DeregnetModel<Model>::addSuboptimalityConstraint(std::set<std::string> nodes_so_far,
                                             double max_overlap,
                                             int size) {
    GRBLinExpr nodes_so_far_expr;
    for (NodeIt v(*graph); v != INVALID; ++v) {
        if (nodes_so_far.find((*nodeid)[v]) != nodes_so_far.end())
            nodes_so_far_expr += x[v];
    }
    model.addConstr(nodes_so_far_expr <= max_overlap * size);
}

template <typename Model>
Solution DeregnetModel<Model>::getCurrentSolution() {
    Solution solution;
    std::set<Node> nodes;
    string rootid;
    for (NodeIt v(*graph); v != INVALID; ++v) {
        if (x[v].get(GRB_DoubleAttr_X) >= 0.98) {
            nodes.insert(v);
            solution.nodes.insert((*nodeid)[v]);
            if (!root) {
                if ((*y)[v].get(GRB_DoubleAttr_X) >= 0.98)
                    rootid = (*nodeid)[v];
            }
            else
                rootid = (*nodeid)[*root];
        }
    }
    for (auto v : nodes) {
        std::string source = (*nodeid)[v];
        for (OutArcIt a(*graph, v); a != INVALID; ++a) {
            std::string target { (*nodeid)[graph->target(a)] };
            if (solution.nodes.find(target) != solution.nodes.end())
                solution.edgelist.insert( std::make_pair(source, target));
        }
    }
    solution.rootid = rootid;
    solution.total_score = model.get(GRB_DoubleAttr_ObjVal);
    solution.avg_score = solution.total_score / solution.nodes.size();
    return solution;
}


// Specializations for Model = GRBModel

template <> inline
void DeregnetModel<GRBModel>::createVariablesRoot() {
    for (NodeIt v(*graph); v != INVALID; ++v)
        x[v] = model.addVar(0.0, 1.0, (*score)[v], GRB_BINARY);
}

template <> inline
void DeregnetModel<GRBModel>::createVariablesNoRoot() {
    y = new std::map<Node, GRBVar>;
    for (NodeIt v(*graph); v != INVALID; ++v) {
        x[v] = model.addVar(0.0, 1.0, (*score)[v], GRB_BINARY);
        (*y)[v] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
    }
}

template <> inline
void DeregnetModel<GRBModel>::setCallbackRoot(double* gap_cut) {
    model.setCallback( new LazyConstraintCallbackRoot(&x, graph, root, gap_cut) );
}

template <> inline
void DeregnetModel<GRBModel>::setCallbackNoRoot(double* gap_cut) {
    model.setCallback( new LazyConstraintCallbackNoRoot(&x, y, graph, gap_cut) );
}

template <> inline
bool DeregnetModel<GRBModel>::solve(std::pair<Node, std::set<Node>>* start_solution,
                                    double* time_limit,
                                    double* gap_cut,
                                    std::string model_sense,
                                    grbfrc::Algorithm algorithm) {
    solve_common(start_solution, time_limit, gap_cut, model_sense);
    model.optimize();
    int status = model.get(GRB_IntAttr_Status);   // FMILP implmentation of status!
    if (status == GRB_OPTIMAL || status == GRB_INTERRUPTED || status == GRB_TIME_LIMIT)
        return true;
    else
        return false;
}

template <> inline
void DeregnetModel<GRBModel>::addSuboptimalityConstraint(std::set<std::string> nodes_so_far,
                                             double max_overlap,
                                             int size) {
    GRBLinExpr nodes_so_far_expr;
    for (NodeIt v(*graph); v != INVALID; ++v) {
        if (nodes_so_far.find((*nodeid)[v]) != nodes_so_far.end())
            nodes_so_far_expr += x[v];
    }
    model.addConstr(nodes_so_far_expr <= max_overlap * size);
}

// Specializations for Model = grbfrc::FMILP

template <> inline
void DeregnetModel<grbfrc::FMILP>::createVariablesRoot() {
    for (NodeIt v(*graph); v != INVALID; ++v)
        x[v] = model.addVar(0.0, 1.0, (*score)[v], 1.0, GRB_BINARY);
}

template <> inline
void DeregnetModel<grbfrc::FMILP>::createVariablesNoRoot() {
    y = new std::map<Node, GRBVar>;
    for (NodeIt v(*graph); v != INVALID; ++v) {
        x[v] = model.addVar(0.0, 1.0, (*score)[v], 1.0, GRB_BINARY);
        (*y)[v] = model.addVar(0.0, 1.0, 0.0, 0.0, GRB_BINARY);
    }
}

template <> inline
void DeregnetModel<grbfrc::FMILP>::setCallbackRoot(double* gap_cut) {

}

template <> inline
void DeregnetModel<grbfrc::FMILP>::setCallbackNoRoot(double* gap_cut) {

}

template <> inline
bool DeregnetModel<grbfrc::FMILP>::solve(std::pair<Node, std::set<Node>>* start_solution,
                                    double* time_limit,
                                    double* gap_cut,
                                    std::string model_sense,
                                    grbfrc::Algorithm algorithm) {
    return false;
    /*
    solve_common(start_solution, time_limit, gap_cut, model_sense);
    model.optimize(algorithm);
    int status = model.get(GRB_IntAttr_Status);   // FMILP implmentation of status!
    if (status == GRB_OPTIMAL || status == GRB_INTERRUPTED || status == GRB_TIME_LIMIT)
        return true;
    else
        return false;
    */
}

template <> inline
void DeregnetModel<grbfrc::FMILP>::addSuboptimalityConstraint(std::set<std::string> nodes_so_far,
                                                              double max_overlap,
                                                              int size) {
    GRBLinExpr nodes_so_far_expr;
    GRBLinExpr size_expr;
    for (NodeIt v(*graph); v != INVALID; ++v) {
        size_expr += x[v];
        if (nodes_so_far.find((*nodeid)[v]) != nodes_so_far.end())
            nodes_so_far_expr += x[v];
    }
    model.addConstr(nodes_so_far_expr <= max_overlap * size_expr);
}

}   // namespace deregnet

#endif    // DEREGNET_MODEL_H
