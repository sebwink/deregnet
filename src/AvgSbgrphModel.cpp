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

#include <iostream>
#include <map>
#include <grbfrc/FMILP.h>
#include <deregnet/usinglemon.h>

#include <deregnet/AvgSbgrphModel.h>

namespace deregnet {

AvgSbgrphModel::AvgSbgrphModel(Graph* graphp,
                               NodeMap<double>* scorep,
                               NodeMap<string>* nodeidp,
                               Node* rootp)
 : graph { graphp },
   score { scorep },
   nodeid { nodeidp },
   root { rootp }
{ }


void AvgSbgrphModel::createVariables() {
    if (root)
        createVariablesRoot();
    else
        createVariablesNoRoot();
    filp.update();
}

void AvgSbgrphModel::createVariablesRoot() {
    for (NodeIt v(*graph); v != INVALID; ++v)
        x[v] = filp.addVar(0.0, 1.0, (*score)[v], 1.0, GRB_BINARY);
}

void AvgSbgrphModel::createVariablesNoRoot() {
    y = new std::map<Node, GRBVar>;
    for (NodeIt v(*graph); v != INVALID; ++v) {
        x[v] = filp.addVar(0.0, 1.0, (*score)[v], 1.0, GRB_BINARY);
        (*y)[v] = filp.addVar(0.0, 1.0, 0.0, 0.0, GRB_BINARY);
    }
}

void AvgSbgrphModel::reverseGraph() {

}

void AvgSbgrphModel::addBaseConstraints(int min_size,
                                        int max_size) {
    if (root)
        addBaseConstraintsRoot(min_size, max_size);
    else
        addBaseConstraintsNoRoot(min_size, max_size);
}

void AvgSbgrphModel::addBaseConstraintsRoot(int min_size,
                                            int max_size) {
    GRBLinExpr subgraph_size_lhs;
    for (NodeIt v(*graph); v != INVALID; ++v) {
        subgraph_size_lhs += x[v];
        GRBLinExpr parent_sum;
        for (InArcIt a(*graph,v); a != INVALID; ++a)
            parent_sum += x[graph->source(a)];
        filp.addConstr(x[v] - parent_sum <= 0);
    }
    filp.addConstr(subgraph_size_lhs <= max_size);
    filp.addConstr(subgraph_size_lhs >= min_size);
    filp.update();
}

void AvgSbgrphModel::addBaseConstraintsNoRoot(int min_size,
                                              int max_size) {
    GRBLinExpr subgraph_size_lhs;
    GRBLinExpr exactly_one_root_lhs;
    for (NodeIt v(*graph); v != INVALID; ++v) {
        filp.addConstr((*y)[v] <= x[v]);
        subgraph_size_lhs += x[v];
        exactly_one_root_lhs += (*y)[v];
        GRBLinExpr parent_sum;
        for (InArcIt a(*graph,v); a != INVALID; ++a)
            parent_sum += x[graph->source(a)];
        filp.addConstr(x[v] - (*y)[v] - parent_sum <= 0);
    }
    filp.addConstr(subgraph_size_lhs <= max_size);
    filp.addConstr(subgraph_size_lhs >= min_size);
    filp.addConstr(exactly_one_root_lhs == 1);
    filp.update();
}

void AvgSbgrphModel::addIncludeConstraints(std::set<Node>* include) {
    for (Node v : *include)
        filp.addConstr(x[v] == 1);
}

void AvgSbgrphModel::addExcludeConstraints(std::set<Node>* exclude) {
    for (Node v : *exclude)
        filp.addConstr(x[v] == 0);
}

void AvgSbgrphModel::addReceptorConstraints(std::set<Node>* receptors) {
    if (!root) {
        for (NodeIt v(*graph); v != INVALID; ++v) {
            if (receptors->find(v) != receptors->end())
                filp.addConstr((*y)[v] <= 1);
            else
                filp.addConstr((*y)[v] == 0);
        }
    }
}

void AvgSbgrphModel::addTerminalConstraints(std::set<Node>* terminals) {

}

bool AvgSbgrphModel::solve(grbfrc::Algorithm algorithm,
                           bool start_heuristic,
                           double* time_limit,
                           double* gap_cut,
                           std::string model_sense) {
    if (model_sense == "min")
        filp.set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);
    else
        filp.set(GRB_IntAttr_ModelSense, GRB_MAXIMIZE);
    std::cout << "Yippie!!!" << std::endl;
    return false;

}

void AvgSbgrphModel::addSuboptimalityConstraint(std::set<std::string> nodes_so_far,
                                                double max_overlap) {

}

AvgSbgrphModel::Solution AvgSbgrphModel::getCurrentSolution() {
    Solution solution;

    return solution;
}

}   // namespace deregnet
