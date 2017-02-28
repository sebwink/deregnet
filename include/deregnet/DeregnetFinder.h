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
#include <functional>

#include <gurobi_c++.h>

#include <grbfrc/FMILP.h>

#include <deregnet/utils.h>
#include <deregnet/usinglemon.h>
#include <deregnet/DeregnetModel.h>
#include <deregnet/StartHeuristic.h>
#include <deregnet/SuboptimalStartHeuristic.h>
#include <deregnet/AvgStartHeuristic.h>
#include <deregnet/AvgSuboptimalStartHeuristic.h>

namespace deregnet {

using FMILP = grbfrc::FMILP;

template<typename ModelType, typename Data>
class DeregnetFinder {

  private:

    Data* data;
    DeregnetModel<ModelType, Data>* model;

  public:

    DeregnetFinder(Data* xdata);
    std::vector<Subgraph> run();

  private:

    void find_start_solution(std::pair<Node, std::set<Node>>** start_solution);
    void find_suboptimal_start_solution(std::pair<Node, std::set<Node>>** start_solution,
                                        std::set<std::string>* nodes_so_far);
    Subgraph to_subgraph(Solution solution,
                         string signature);
    void run_optimal_init(std::pair<Node, std::set<Node>>** start_solution);
    void run_suboptimal_init(std::pair<Node, std::set<Node> > **start_solution,
                             std::set<std::string>* nodes_so_far);
    void run_optimal_windup(bool solve_successful,
                            std::vector<Subgraph>* subgraphs);
    bool run_suboptimal_windup(bool solve_successful,
                               std::vector<Subgraph>* subgraphs,
                               std::set<std::string>* nodes_so_far,
                               int i);

};

template <typename ModelType, typename Data>
DeregnetFinder<ModelType, Data>::DeregnetFinder(Data* xdata)
 : data { xdata },
   model { new DeregnetModel<ModelType, Data>(xdata) }
{ }

template <typename ModelType, typename Data>
void DeregnetFinder<ModelType, Data>::run_optimal_init(std::pair<Node, std::set<Node> >** start_solution) {
    model->createVariables();
    model->addBaseConstraints();
    if (data->include)
        model->addIncludeConstraints();
    if (data->exclude)
        model->addExcludeConstraints();
    if (data->receptors && !data->root)
        model->addReceptorConstraints();
    if (data->terminals)
        model->addTerminalConstraints();
    if (data->start_heuristic)
        find_start_solution(start_solution);
    if (data->start_heuristic && !(*start_solution))
        std::cout << "No heuristic start solution found.\n" << std::endl;
}

template <typename ModelType, typename Data>
void DeregnetFinder<ModelType, Data>::run_optimal_windup(bool solve_successful, std::vector<Subgraph> *subgraphs) {
    if ( solve_successful ) {
        subgraphs->push_back ( to_subgraph( model->getCurrentSolution(), "optimal" ) );
    }
    else {
        std::cerr << "No optimal subgraph could be found." << std::endl;
        exit(NO_OPTIMAL_SUBGRAPH_FOUND);
    }
}

template <typename ModelType, typename Data>
void DeregnetFinder<ModelType, Data>::run_suboptimal_init(std::pair<Node, std::set<Node> > **start_solution,
                                                          std::set<std::string>* nodes_so_far) {
    model->addSuboptimalityConstraint(*nodes_so_far);
    if (data->start_heuristic)
        find_suboptimal_start_solution(start_solution, nodes_so_far);
    if (data->start_heuristic && !(*start_solution))
        std::cout << "No heuristic start solution found.\n" << std::endl;
}

template <typename ModelType, typename Data>
bool DeregnetFinder<ModelType, Data>::run_suboptimal_windup(bool solve_successful,
                                                            std::vector<Subgraph>* subgraphs,
                                                            std::set<std::string>* nodes_so_far,
                                                            int i) {
    if ( solve_successful ) {
        Solution current = model->getCurrentSolution() ;
        subgraphs->push_back(to_subgraph( current, "suboptimal_" + std::to_string(i) ));
        for (auto node: current.nodes)
            nodes_so_far->insert(node);
        return true;
    }
    else {
        std::cout << "With these parameters, no more than " << i
             << " suboptimal subgraphs could be found." << std::endl;
        return false;
    }
}

template <typename ModelType, typename Data>
std::vector<Subgraph> DeregnetFinder<ModelType, Data>::run() {
    std::pair<Node, std::set<Node>>* start_solution { nullptr };
    // find optimal subgraph
    run_optimal_init(&start_solution);

    /*
    std::cout << "# nodes in start: " << (start_solution->second).size() << std::endl;
    write2sif(data->graph, start_solution->second, data->nodeid, "start.sif");
    */

    std::vector<Subgraph> subgraphs;
    bool solve_successful { model->solve(start_solution) };
    run_optimal_windup(solve_successful, &subgraphs);
    // find suboptimal subgraphs
    std::set<std::string> nodes_so_far { subgraphs.back().nodes };
    for (int i = 0; i < data->num_subopt_iter; ++i) {
        start_solution = nullptr;
        run_suboptimal_init(&start_solution, &nodes_so_far);
        // solve for next suboptimal solution
        solve_successful = model->solve(start_solution);
        if (!run_suboptimal_windup(solve_successful, &subgraphs, &nodes_so_far, i))
            break;
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

// specializations

template <> inline
void DeregnetFinder<GRBModel, DrgntData>::find_start_solution(std::pair<Node,
                                                              std::set<Node>>** start_solution) {
    std::function<bool(double, double)> cmp;
    if (data->model_sense == "min") cmp = std::greater<double>();
    else cmp = std::less<double>();

    StartHeuristic heuristic(data->graph, data->score, data->root, data->exclude, data->receptors, cmp, data->size);
    if (heuristic.run())
        *start_solution = heuristic.getStartSolution();
}


template <> inline
void DeregnetFinder<FMILP, AvgdrgntData>::find_start_solution(std::pair<Node,
                                                              std::set<Node>>** start_solution) {


    std::function<bool(double, double)> cmp;
    if (data->model_sense == "min") cmp = std::greater<double>();
    else cmp = std::less<double>();
    AvgStartHeuristic heuristic(data->graph, data->score, data->root, data->exclude, data->receptors,
                                cmp, data->min_size, data->max_size);
    if (heuristic.run())
        *start_solution = heuristic.getStartSolution();
}


template <> inline
void DeregnetFinder<GRBModel, DrgntData>::find_suboptimal_start_solution(std::pair<Node, std::set<Node>>** start_solution,
                                                                         std::set<std::string>* nodes_so_far) {
    std::function<bool(double, double)> cmp;
    if (data->model_sense == "min") cmp = std::greater<double>();
    else cmp = std::less<double>();
    SuboptimalStartHeuristic heuristic(data->graph, data->score, data->root,
                                       data->exclude, data->receptors, cmp,
                                       data->nodeid, nodes_so_far, data->max_overlap, data->size);
    if (heuristic.run())
        *start_solution = heuristic.getStartSolution();
}


template <> inline
void DeregnetFinder<FMILP, AvgdrgntData>::find_suboptimal_start_solution(std::pair<Node, std::set<Node>>** start_solution,
                                                                         std::set<std::string>* nodes_so_far) {

    std::function<bool(double, double)> cmp;
    if (data->model_sense == "min") cmp = std::greater<double>();
    else cmp = std::less<double>();
    AvgSuboptimalStartHeuristic heuristic(data->graph, data->score, data->root,
                                          data->exclude, data->receptors, cmp,
                                          data->nodeid, nodes_so_far, data->max_overlap,
                                          data->min_size, data->max_size);
    if (heuristic.run())
        *start_solution = heuristic.getStartSolution();

}


}

#endif    //    DEREGNET_FINDER_H
