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
#include <deregnet/GrbfrcLazyConstraintCallback.h>
#include <deregnet/DrgntData.h>
#include <deregnet/AvgdrgntData.h>

namespace deregnet {

template <typename Model, typename Data>
class DeregnetModel {

  protected:

    GRBEnv env;
    Model model { Model(env) };
    std::map<Node, GRBVar> x;
    std::map<Node, GRBVar>* y { nullptr };

    Data* data;
    Graph* graph;
    NodeMap<double>* score;
    NodeMap<std::string>* nodeid;
    Node* root;


  public:

    DeregnetModel(Data* data);
    void createVariables();
    void addBaseConstraints();
    void addIncludeConstraints();
    void addExcludeConstraints();
    void addReceptorConstraints();
    void addTerminalConstraints();
    bool solve(std::pair<Node, std::set<Node>>* start_solution);
    void addSuboptimalityConstraint(std::set<std::string>& nodes_so_far);
    Solution getCurrentSolution();
    void setObjSense(int sense);

  private:

    void createVariablesRoot();
    void createVariablesNoRoot();
    void addBaseConstraintsRoot();
    void addBaseConstraintsNoRoot();
    void setStartSolution(std::pair<Node, std::set<Node>>* start_solution);  // specialize for GRBModel and FMILP ?
    void setup_solve(std::pair<Node, std::set<Node>>* start_solution);
    void setCallbackRoot();
    void setCallbackNoRoot();
    void _getCurrentSolution(std::string* rootid, std::set<Node>* nodes, Solution* solution);

    GRBLinExpr setBaseConstraintsRootCommon();
    GRBLinExpr setBaseConstraintsNoRootCommon();

};

template <typename Model, typename Data>
DeregnetModel<Model, Data>::DeregnetModel(Data* xdata)
 : data { xdata },
   graph { xdata->graph },
   score { xdata->score },
   nodeid { xdata->nodeid },
   root { xdata->root }
{ }

template <typename Model, typename Data>
void DeregnetModel<Model, Data>::createVariables() {
    if (root)
        createVariablesRoot();
    else
        createVariablesNoRoot();
    model.update();
}

template <typename Model, typename Data>
void DeregnetModel<Model, Data>::addBaseConstraints() {
    if (root)
        addBaseConstraintsRoot();
    else
        addBaseConstraintsNoRoot();
    model.update();
}

template <typename Model, typename Data>
GRBLinExpr DeregnetModel<Model, Data>::setBaseConstraintsRootCommon() {
    GRBLinExpr subgraph_size_lhs;
    for (NodeIt v(*graph); v != INVALID; ++v) {
        subgraph_size_lhs += x[v];
        GRBLinExpr parent_sum;
        for (InArcIt a(*graph,v); a != INVALID; ++a)
            parent_sum += x[graph->source(a)];
        if (v != *root)
            model.addConstr(x[v] - parent_sum <= 0);
    }
    model.addConstr(x[*root] == 1);
    return subgraph_size_lhs;
}

template <typename Model, typename Data>
GRBLinExpr DeregnetModel<Model, Data>::setBaseConstraintsNoRootCommon() {
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
    model.addConstr(exactly_one_root_lhs == 1);
    return subgraph_size_lhs;
}

template <> inline
void DeregnetModel<GRBModel, DrgntData>::addBaseConstraintsRoot() {
    GRBLinExpr subgraph_size_lhs { setBaseConstraintsRootCommon() };
    model.addConstr(subgraph_size_lhs == data->size);
}

template <> inline
void DeregnetModel<GRBModel, DrgntData>::addBaseConstraintsNoRoot() {
    GRBLinExpr subgraph_size_lhs { setBaseConstraintsNoRootCommon() };
    model.addConstr(subgraph_size_lhs == data->size);
}

template <> inline
void DeregnetModel<grbfrc::FMILP, AvgdrgntData>::addBaseConstraintsRoot() {
    GRBLinExpr subgraph_size_lhs { setBaseConstraintsRootCommon() };
    model.addConstr(data->min_size <= subgraph_size_lhs);
    model.addConstr(subgraph_size_lhs <= data->max_size);
}

template <> inline
void DeregnetModel<grbfrc::FMILP, AvgdrgntData>::addBaseConstraintsNoRoot() {
    GRBLinExpr subgraph_size_lhs { setBaseConstraintsNoRootCommon() };
    model.addConstr(data->min_size <= subgraph_size_lhs);
    model.addConstr(subgraph_size_lhs <= data->max_size);
}

template <typename Model, typename Data>
void DeregnetModel<Model, Data>::addIncludeConstraints() {
    for (Node v : *(data->include))
        model.addConstr(x[v] == 1);
    model.update();
}

template <typename Model, typename Data>
void DeregnetModel<Model, Data>::addExcludeConstraints() {
    for (Node v : *(data->exclude))
        model.addConstr(x[v] == 0);
    model.update();
}

template <typename Model, typename Data>
void DeregnetModel<Model, Data>::addReceptorConstraints() {
    if (!root) {
        std::set<Node>* receptors { data->receptors };
        for (NodeIt v(*graph); v != INVALID; ++v)
            if (receptors->find(v) == receptors->end()) {
                model.addConstr((*y)[v] == 0);
            }
    }
    model.update();
}

template <> inline
void DeregnetModel<grbfrc::FMILP, AvgdrgntData>::addReceptorConstraints() {
    if (!root) {
        std::set<Node>* receptors { data->receptors };
        GRBLinExpr receptor_sum;
        for (NodeIt v(*graph); v != INVALID; ++v)
            if (receptors->find(v) == receptors->end()) {
                model.addConstr((*y)[v] == 0);   
            }
            else {
                receptor_sum += x[v];
            }
        if (data->min_num_receptors)
            model.addConstr(receptor_sum >= *(data->min_num_receptors));
    }
    model.update();
}

template <typename Model, typename Data>
void DeregnetModel<Model, Data>::addTerminalConstraints() { // consider using complement of terminals ...
    for (NodeIt v(*graph); v != INVALID; ++v) {
        std::set<Node>* terminals { data->terminals };
        if (terminals->find(v) == terminals->end()) {
            GRBLinExpr out_neighbor_expr;
            for (OutArcIt a(*graph, v); a != INVALID; ++a)
                out_neighbor_expr += x[graph->target(a)];
            model.addConstr(x[v] - out_neighbor_expr <= 0);
        }
    }
    model.update();
}

template <> inline
void DeregnetModel<grbfrc::FMILP, AvgdrgntData>::addTerminalConstraints() { // consider using complement of terminals ...
    for (NodeIt v(*graph); v != INVALID; ++v) {
        std::set<Node>* terminals { data->terminals };
        GRBLinExpr terminal_sum;
        if (terminals->find(v) == terminals->end()) {
            terminal_sum += x[v];
            GRBLinExpr out_neighbor_expr;
            for (OutArcIt a(*graph, v); a != INVALID; ++a)
                out_neighbor_expr += x[graph->target(a)];
            model.addConstr(x[v] - out_neighbor_expr <= 0);
        }
        if (data->min_num_terminals)
            model.addConstr(terminal_sum >= *(data->min_num_terminals));
    }
    model.update();
}

template<typename Model, typename Data>
void DeregnetModel<Model, Data>::setObjSense(int sense) {
    model.set(GRB_IntAttr_ModelSense, sense);
}

template <> inline
void DeregnetModel<grbfrc::FMILP, AvgdrgntData>::setObjSense(int sense) {
    model.setObjSense(sense);
}

template <typename Model, typename Data>
void DeregnetModel<Model, Data>::setup_solve(std::pair<Node, std::set<Node>>* start_solution) {
    if (data->model_sense == "min")
        setObjSense(GRB_MINIMIZE);
    else
        setObjSense(GRB_MAXIMIZE);
    model.update();

    model.set(GRB_IntParam_LazyConstraints, 1);

    if (data->time_limit)
        model.set(GRB_DoubleParam_TimeLimit, *(data->time_limit));

    if (start_solution)
        setStartSolution(start_solution);

    if (root)
        setCallbackRoot();
    else
        setCallbackNoRoot();
    model.update();
}


template <typename Model, typename Data>
Solution DeregnetModel<Model, Data>::getCurrentSolution() {
    Solution solution;
    std::set<Node> nodes;
    string rootid;
    _getCurrentSolution(&rootid, &nodes, &solution);
    for (auto v : nodes) {
        std::string source = (*nodeid)[v];
        for (OutArcIt a(*graph, v); a != INVALID; ++a) {
            std::string target { (*nodeid)[graph->target(a)] };
            if (solution.nodes.find(target) != solution.nodes.end())
                solution.edgelist.insert( std::make_pair(source, target));
        }
    }
    solution.rootid = rootid;
    solution.avg_score = solution.total_score / solution.nodes.size();
    return solution;
}


// Specializations

template <> inline
void DeregnetModel<GRBModel, DrgntData>::createVariablesRoot() {
    for (NodeIt v(*graph); v != INVALID; ++v)
        x[v] = model.addVar(0.0, 1.0, (*score)[v], GRB_BINARY);
}

template <> inline
void DeregnetModel<grbfrc::FMILP, AvgdrgntData>::createVariablesRoot() {
    for (NodeIt v(*graph); v != INVALID; ++v)
        x[v] = model.addVar(0.0, 1.0, (*score)[v], 1.0, GRB_BINARY);
}

template <> inline
void DeregnetModel<GRBModel, DrgntData>::createVariablesNoRoot() {
    y = new std::map<Node, GRBVar>;
    for (NodeIt v(*graph); v != INVALID; ++v) {
        x[v] = model.addVar(0.0, 1.0, (*score)[v], GRB_BINARY);
        (*y)[v] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
    }
}

template <> inline
void DeregnetModel<grbfrc::FMILP, AvgdrgntData>::createVariablesNoRoot() {
    y = new std::map<Node, GRBVar>;
    for (NodeIt v(*graph); v != INVALID; ++v) {
        x[v] = model.addVar(0.0, 1.0, (*score)[v], 1.0, GRB_BINARY);
        (*y)[v] = model.addVar(0.0, 1.0, 0.0, 0.0, GRB_BINARY);
    }
}

template <> inline
void DeregnetModel<GRBModel, DrgntData>::addSuboptimalityConstraint(std::set<std::string>& nodes_so_far) {
    GRBLinExpr nodes_so_far_expr;
    for (NodeIt v(*graph); v != INVALID; ++v) {
        if (nodes_so_far.find((*nodeid)[v]) != nodes_so_far.end())
            nodes_so_far_expr += x[v];
    }
    model.addConstr(nodes_so_far_expr <= data->max_overlap * data->size);
}

template <> inline
void DeregnetModel<grbfrc::FMILP, AvgdrgntData>::addSuboptimalityConstraint(std::set<std::string>& nodes_so_far) {
    GRBLinExpr nodes_so_far_expr;
    GRBLinExpr size_expr;
    for (NodeIt v(*graph); v != INVALID; ++v) {
        size_expr += x[v];
        if (nodes_so_far.find((*nodeid)[v]) != nodes_so_far.end())
            nodes_so_far_expr += x[v];
    }
    model.addConstr(nodes_so_far_expr <= data->max_overlap * size_expr);
}

template <> inline
void DeregnetModel<GRBModel, DrgntData>::setStartSolution(std::pair<Node, std::set<Node>>* start_solution) {
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

template <> inline
void DeregnetModel<grbfrc::FMILP, AvgdrgntData>::setStartSolution(std::pair<Node, std::set<Node>>* start_solution) {
    for (NodeIt v(*graph); v != INVALID; ++v) {
        if ((start_solution->second).find(v) != (start_solution->second).end())
            model.setStartSolution(x[v], 1.0);
        else
            model.setStartSolution(x[v], 1.0);
        if (!root) {
            if (v == start_solution->first)
                model.setStartSolution((*y)[v], 1.0);
            else
                model.setStartSolution((*y)[v], 0.0);
        }
    }
}


template <> inline
void DeregnetModel<GRBModel, DrgntData>::setCallbackRoot() {
    model.setCallback( new LazyConstraintCallbackRoot(&x, graph, root, data->gap_cut) );
}

template <> inline
void DeregnetModel<grbfrc::FMILP, AvgdrgntData>::setCallbackRoot() {
}

template <> inline
void DeregnetModel<GRBModel, DrgntData>::setCallbackNoRoot() {
    model.setCallback( new LazyConstraintCallbackNoRoot(&x, y, graph, data->gap_cut) );
}

template <> inline
void DeregnetModel<grbfrc::FMILP, AvgdrgntData>::setCallbackNoRoot() {
}

template <> inline
bool DeregnetModel<GRBModel, DrgntData>::solve(std::pair<Node, std::set<Node>>* start_solution) {
    setup_solve(start_solution);
    model.optimize();
    int status = model.get(GRB_IntAttr_Status);   // FMILP implmentation of status!
    return (status == GRB_OPTIMAL || status == GRB_INTERRUPTED || status == GRB_TIME_LIMIT);
}

template <> inline
bool DeregnetModel<grbfrc::FMILP, AvgdrgntData>::solve(std::pair<Node, std::set<Node>>* start_solution) {
    setup_solve(start_solution);
    // add callback here
    grbfrc::GrbfrcCallback<Node>* cb;
    if (root) {
        std::vector<std::map<Node, GRBVar>>* vargrps = new std::vector<std::map<Node, GRBVar>>( { x } );
        GrbfrcLazyConstraintCallbackRoot* callback = new GrbfrcLazyConstraintCallbackRoot(graph, root, data->gap_cut);
        cb = new grbfrc::GrbfrcCallback<Node>(callback, vargrps);
    }
    else {
        std::vector<std::map<Node, GRBVar>>* vargrps = new std::vector<std::map<Node, GRBVar>>( { x, *y } );
        GrbfrcLazyConstraintCallbackNoRoot* callback = new GrbfrcLazyConstraintCallbackNoRoot(graph, data->gap_cut);
        cb = new grbfrc::GrbfrcCallback<Node>(callback, vargrps);
    }
    double objlb = 0.0;
    double objub = 0.0;
    std::set<Node> maxnodes;
    for (int k = 0; k <= data->min_size; ++k) {
        double ub = -1000000.0;
        for (NodeIt v(*graph); v != INVALID; ++v) {
            if (maxnodes.find(v) == maxnodes.end() && score->operator[](v) > ub) {
                ub = score->operator[](v);
                objub += ub;
                maxnodes.insert(v);
            }
        }
    }

    objub = objub / data->min_size;
    // std::cout << objub << std::endl;

    model.optimize(data->algorithm, cb, &objub, &objlb);
    int status = model.get(GRB_IntAttr_Status);
    return (status == GRB_OPTIMAL || status == GRB_INTERRUPTED || status == GRB_TIME_LIMIT);
}

template <> inline
void DeregnetModel<GRBModel, DrgntData>::_getCurrentSolution(std::string* rootid, std::set<Node>* nodes, Solution* solution) {
    for (NodeIt v(*graph); v != INVALID; ++v) {
        if (x[v].get(GRB_DoubleAttr_X) >= 0.98) {
            nodes->insert(v);
            (solution->nodes).insert((*nodeid)[v]);
            if (!root) {
                if ((*y)[v].get(GRB_DoubleAttr_X) >= 0.98)
                    *rootid = (*nodeid)[v];
            }
            else
                *rootid = (*nodeid)[*root];
        }
    }
    solution->total_score = model.get(GRB_DoubleAttr_ObjVal);
}

template <> inline
void DeregnetModel<grbfrc::FMILP, AvgdrgntData>::_getCurrentSolution(std::string* rootid, std::set<Node>* nodes, Solution* solution) {
    for (NodeIt v(*graph); v != INVALID; ++v) {
        if (model.getVal(x[v]) >= 0.98) {
            nodes->insert(v);
            (solution->nodes).insert((*nodeid)[v]);
            if (!root) {
                if (model.getVal((*y)[v]) >= 0.98)
                    *rootid = (*nodeid)[v];
            }
            else
                *rootid = (*nodeid)[*root];
        }
    }
    solution->total_score = model.getObjVal() * nodes->size();
}


}   // namespace deregnet

#endif    // DEREGNET_MODEL_H
