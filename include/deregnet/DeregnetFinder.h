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

#ifndef DEREGNET_FINDER_H
#define DEREGNET_FINDER_H

#include <iostream>
#include <string>
#include <vector>

#include <gurobi_c++.h>

#include <deregnet/utils.h>
#include <deregnet/usinglemon.h>
#include <deregnet/DrgntData.h>
#include <deregnet/DeregnetModel.h>
#include <deregnet/StartHeuristic.h>
#include <deregnet/SuboptimalStartHeuristic.h>

namespace deregnet {

template<typename ModelType, typename Data>
class DeregnetFinder {

  private:

    Data* data;
    DeregnetModel<ModelType>* model;

  public:

    DeregnetFinder(Data* xdata);
    std::vector<Subgraph> run(bool start_heuristic, std::string model_sense);

  private:

    void find_start_solution(std::pair<Node, std::set<Node>>** start_solution, std::string model_sense);
    void find_suboptimal_start_solution(std::pair<Node, std::set<Node>>** start_solution,
                                        std::string model_sense,
                                        std::set<std::string>* nodes_so_far);
    void add_base_constraints();
    void add_suboptimiality_constraints(std::set<std::string>& nodes_so_far);
    Subgraph to_subgraph(Solution solution, string signature);

};

template <typename ModelType, typename Data>
DeregnetFinder<ModelType, Data>::DeregnetFinder(Data* xdata)
 : data { xdata },
   model { new DeregnetModel<ModelType>(xdata->graph, xdata->score, xdata->nodeid, xdata->root) }
{ }

template <typename ModelType, typename Data>
std::vector<Subgraph> DeregnetFinder<ModelType, Data>::run(bool start_heuristic, std::string model_sense) {
    model->createVariables();
    add_base_constraints();
    if (data->include)
        model->addIncludeConstraints(data->include);
    if (data->exclude)
        model->addExcludeConstraints(data->exclude);
    if (data->receptors && !data->root)
        model->addReceptorConstraints(data->receptors);
    if (data->terminals)
        model->addTerminalConstraints(data->terminals);

    std::vector<Subgraph> subgraphs;
    // find optimal subgraph (modulo time_limit and/or gap_cut)

    std::pair<Node, std::set<Node>>* start_solution { nullptr };
    if (start_heuristic)
        find_start_solution(&start_solution, model_sense);

    if (!start_solution)
        std::cout << "No heuristic start solution found.\n" << std::endl;

    if ( model->solve(start_solution,
                     data->time_limit,
                     data->gap_cut,
                     model_sense) ) {
        subgraphs.push_back ( to_subgraph( model->getCurrentSolution(), "optimal" ) );
    }
    else {
        std::cerr << "No optimal subgraph could be found." << std::endl;
        exit(NO_OPTIMAL_SUBGRAPH_FOUND);
    }

    // find suboptimal subgraphs
    std::set<std::string> nodes_so_far { subgraphs.back().nodes };
    for (int i = 0; i < data->num_subopt_iter; i++) {
        add_suboptimiality_constraints(nodes_so_far);
        // suboptimal start solution
        start_solution = nullptr;
        if (start_heuristic)
            find_suboptimal_start_solution(&start_solution, model_sense, &nodes_so_far);

        if (start_heuristic && !start_solution)
            std::cout << "No heuristic start solution found.\n" << std::endl;

        // solve for next suboptimal solution
        if ( model->solve(start_solution,
                         data->time_limit,
                         data->gap_cut,
                         model_sense) ) {
            Solution current = model->getCurrentSolution() ;
            subgraphs.push_back(to_subgraph( current, "suboptimal_" + std::to_string(i) ));
            for (auto node: current.nodes)
                nodes_so_far.insert(node);
        }
        else {
            std::cout << "With these parameters, no more than " << i
                 << " suboptimal subgraphs could be found." << std::endl;
            break;
        }
    }

    return subgraphs;
}

template <typename ModelType, typename Data>
Subgraph DeregnetFinder<ModelType, Data>::to_subgraph(Solution solution, std::string signature) {
    Subgraph subgraph;
    subgraph.signature = signature;
    subgraph.rootid = solution.rootid;
    subgraph.nodes = solution.nodes;
    subgraph.total_score = solution.total_score;
    subgraph.avg_score = solution.avg_score;
    if (!data->receptor_as_root)
        swap(data->receptors, data->terminals);
    if (data->receptors)
        for (auto nodeid : subgraph.nodes) {
            Node node;
            if (getNodeById(data->graph, data->nodeid, nodeid, &node))
                if (data->receptors->find(node) != data->receptors->end())
                    subgraph.receptors.insert(nodeid);
        }
    if (data->terminals)
        for (auto nodeid : subgraph.nodes) {
            Node node;
            if (getNodeById(data->graph, data->nodeid, nodeid, &node))
                if (data->terminals->find(node) != data->terminals->end())
                    subgraph.terminals.insert(nodeid);
        }
    if (!data->receptor_as_root)
        for (auto edge : solution.edgelist)
            subgraph.edgelist.insert(std::make_pair(edge.second, edge.first));
    else
        subgraph.edgelist = solution.edgelist;
    return subgraph;
}

// ModelType = GRBModel specializations

template <> inline
void DeregnetFinder<GRBModel, DrgntData>::find_start_solution(std::pair<Node,
                                                              std::set<Node>>** start_solution,
                                                              std::string model_sense) {
    StartHeuristic* heuristic;
    if (model_sense == "min")
        heuristic = new StartHeuristic(data->graph, data->score, data->root, data->size, data->exclude, data->receptors, std::greater<double>());
    else
        heuristic = new StartHeuristic(data->graph, data->score, data->root, data->size, data->exclude, data->receptors, std::less<double>());
    if (heuristic->run())
        *start_solution = heuristic->getStartSolution();
}

template <> inline
void DeregnetFinder<GRBModel, DrgntData>::find_suboptimal_start_solution(std::pair<Node, std::set<Node>>** start_solution,
                                                                         std::string model_sense,
                                                                         std::set<std::string>* nodes_so_far) {
    SuboptimalStartHeuristic* heuristic;
    if (model_sense == "min")
        heuristic = new SuboptimalStartHeuristic(data->graph, data->score, data->root, data->size,
                                                 data->exclude, data->receptors, std::greater<double>(),
                                                 data->nodeid, nodes_so_far, data->max_overlap);
    else
        heuristic = new SuboptimalStartHeuristic(data->graph, data->score, data->root, data->size,
                                                 data->exclude, data->receptors, std::less<double>(),
                                                 data->nodeid, nodes_so_far, data->max_overlap);
    if (heuristic->run())
        *start_solution = heuristic->getStartSolution();
}

template <> inline
void DeregnetFinder<GRBModel, DrgntData>::add_base_constraints() {
    model->addBaseConstraints(data->size);
}

template <> inline
void DeregnetFinder<GRBModel, DrgntData>::add_suboptimiality_constraints(std::set<std::string>& nodes_so_far) {
    model->addSuboptimalityConstraint(nodes_so_far, data->max_overlap, data->size);
}

}

#endif    //    DEREGNET_FINDER_H
