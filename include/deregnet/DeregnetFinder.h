// --------------------------------------------------------------------------
//                   deregnet -- find deregulated pathways
// --------------------------------------------------------------------------
// Copyright Sebastian Winkler --- Eberhard Karls University Tuebingen, 2018
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

/** 
 * @file DeregnetFinder.h
 *
 * @brief Declares the DeregnetFinder class.
 *
 * @author Sebastian Winkler
 *
 * The DeregnetFinder class provides the main entry point
 * for using the DeRegNet algorithms in C++.
 */

#ifndef DEREGNET_FINDER_H
#define DEREGNET_FINDER_H

#include <fstream>

#include <deregnet/DeregnetModel.h>
#include <deregnet/StartHeuristic.h>
#include <deregnet/SuboptimalStartHeuristic.h>
#include <deregnet/AvgStartHeuristic.h>
#include <deregnet/AvgSuboptimalStartHeuristic.h>

namespace deregnet {

using FMILP = grbfrc::FMILP;

/**
 * @brief Main entry point for running DeRegNet algorithms.
 *
 * DeregnetFinder has two template paramters: 1) ModelType,
 * referring to either GRBModel or grbfrc::FMILP, 2) Data, 
 * referring to a type representing algorithm data suitable
 * for the ModelType parameter. In case of ModelType = GRBModel,
 * Data should be DrgntData, in case of ModelType = grbfrc::FMILP,
 * Data should be AvgdrgntData.
 *
 * @todo This desgin seems not very optimal and somehow flawed. The 
 * most straight-forward improvment would be to implement the
 * absolute version of the algorithm, i.e. ModelType = GRBModel, in
 * terms of grbfrc::FMILP. ModelType would go away and Data could specfiy
 * which algorithm to run, for example.
 *
 * @todo More refined logging system.
 *
 * @author Sebastian Winkler
 */
template<typename ModelType, typename Data>
class DeregnetFinder {

  private:

    Data* data; ///< Pointer to DeregnetData, i.e. DrgntData or AvgdrgntData
    DeregnetModel<ModelType, Data>* model; ///< Model instance to be build and solved

  public:

    /**
     * @brief Constructor, taking the relevant data as argument
     *
     * A DeregnetFinder object is instantiated by passing a pointer to
     * a suitable Data object. In case of ModelType = GRBModel, this
     * data object should be of type DrgntData, in case of ModelType =
     * grbfrc::FMILP, it should be of type AvgdrngtData. (Having to
     * state this here is most likely expression of a design flaw.
     * As long as you do as just specfified, it works though and otherwise
     * the code will presumably not compile anyhow. But having to match the
     * two template parameters manually feels wrong.)
     *
     * @param[in] xdata Pointer to suitable DeregnetData object. If ModelType = GRBModel,
     * you have to pass a pointer to a DrgntData object, in case of ModelType =
     * grbfrc::FMILP, you have to pass a pointer to a AvgdrgntData object.
     */
    DeregnetFinder(Data* xdata);

    /**
     * @brief Method which finds the subgraphs.
     *
     * This method finds subgraphs as specified on instantiation of the DeregnetFinder
     * instance. You get back a vector with the subgraphs.
     *
     * @return Vector containing the the subgraphs
     */
    std::vector<Subgraph> run();

    /**
     * @brief Finds a heuristic optimal start solution.
     */
    void find_start_solution(std::pair<Node, std::set<Node>>** start_solution);
    
    /**
     * @brief Finds a heuristic suboptimal start solution.
     */
    void find_suboptimal_start_solution(std::pair<Node, std::set<Node>>** start_solution,
                                        std::set<std::string>* nodes_so_far);

    /**
     * @brief Transform a solution from the Solution type to the Subgraph type.
     */
    Subgraph to_subgraph(Solution solution,
                         std::string signature);

  private:

    /**
     * @brief Preparation for finding the optimal subgraph
     */
    void run_optimal_init(std::pair<Node, std::set<Node>>** start_solution);

    /**
     * @brief Perparation for finding suboptimal subgraphs
     */
    void run_suboptimal_init(std::pair<Node, std::set<Node> > **start_solution,
                             std::set<std::string>* nodes_so_far);

    /**
     * @brief Post-processing after attempt to find optimal subgraph
     */
    void run_optimal_windup(bool solve_successful,
                            std::vector<Subgraph>* subgraphs);

    /**
     * @brief Post-processing after attempt to find suboptimal subgraphs
     */
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
    model->createVariables();       // create binary variables corresponding to the nodes
    model->addBaseConstraints();    // add base constraints to the model, ie connectivity etc.
    // add optional constraints (ie. include, exclude, receptor and/or terminal constraints)
    if (data->include)
        model->addIncludeConstraints();   // which nodes to include in the subgraph, no matter what
    if (data->exclude)
        model->addExcludeConstraints();   // which nodes to exclude from the subgraph, no matter what
    if (data->receptors && !data->root)
        model->addReceptorConstraints();  // how to handle nodes which are 'receptors'
    if (data->terminals)
        model->addTerminalConstraints();  // how to handle nodes which are 'temrinals'
    // optionally find heuristic start solution by greedy algorithm
    if (data->start_heuristic)
        find_start_solution(start_solution);
    // in case a heuristic start solution was specfied but none was found, issue a warning.
    // TODO: not necessarily always to std::cout ...
    if (data->start_heuristic && !(*start_solution))
        std::cout << "No heuristic start solution found.\n" << std::endl;
}


template <typename ModelType, typename Data>
void DeregnetFinder<ModelType, Data>::run_optimal_windup(bool solve_successful, std::vector<Subgraph> *subgraphs) {
    // if a solution was found update subgraphs data structure
    if ( solve_successful ) {
        subgraphs->push_back ( to_subgraph( model->getCurrentSolution(), "optimal" ) );
    }
    // if no solution was found, issue a message and end the program
    // TODO: not necessarily always to std::cerr
    // TODO: exit normally and handle issue in normal program flow?
    else {
        std::cerr << "No optimal subgraph could be found." << std::endl;
        exit(NO_OPTIMAL_SUBGRAPH_FOUND);
    }
}


template <typename ModelType, typename Data>
void DeregnetFinder<ModelType, Data>::run_suboptimal_init(std::pair<Node, std::set<Node> > **start_solution,
                                                          std::set<std::string>* nodes_so_far) {
    // make sure only a certain fraction of nodes contained in some already found subgraph
    // will be contained this subgraph
    model->addSuboptimalityConstraint(*nodes_so_far);  
    // optionally find heuristic start solution by greedy algorithm
    if (data->start_heuristic)
        find_suboptimal_start_solution(start_solution, nodes_so_far);
    // in case a heuristic start solution was specfied but none was found, issue a warning.
    // TODO: not necessarily always to std::cout ...
    if (data->start_heuristic && !(*start_solution))
        std::cout << "No heuristic start solution found.\n" << std::endl;
}


template <typename ModelType, typename Data>
bool DeregnetFinder<ModelType, Data>::run_suboptimal_windup(bool solve_successful,
                                                            std::vector<Subgraph>* subgraphs,
                                                            std::set<std::string>* nodes_so_far,
                                                            int i) {
    // if a solution was found update subgraphs data structure
    // Also, update the data structure tracking nodes contained in any already found subgraph
    if ( solve_successful ) {
        Solution current = model->getCurrentSolution() ;
        subgraphs->push_back(to_subgraph( current, "suboptimal_" + std::to_string(i) ));
        for (auto node: current.nodes)
            nodes_so_far->insert(node);
        return true;
    }
    // if no solution was found, issue a message
    // TODO: not necessarily always to std::cout
    else {
        std::cout << "With these parameters, no more than " << i
                  << " suboptimal subgraphs could be found." << std::endl;
        return false;
    }
}


template <typename ModelType, typename Data>
std::vector<Subgraph> DeregnetFinder<ModelType, Data>::run() {
    std::pair<Node, std::set<Node>>* start_solution { nullptr };      // potential heuristic start solution
    // find optimal subgraph
    run_optimal_init(&start_solution);                                // construct model, etc.
/*
    std::ofstream file;
    file.open("heuristic_start.nodes");
    for (auto v : start_solution->second)
        file << (*data->nodeid)[v] << "\n";
    file.close();
*/
    std::vector<Subgraph> subgraphs;                                  // data structure for found subgraphs
    bool solve_successful { model->solve(start_solution) };           // try to find optimal solution
    run_optimal_windup(solve_successful, &subgraphs);                 // handle result of optimal solution attempt
    // find suboptimal subgraphs
    std::set<std::string> nodes_so_far { subgraphs.back().nodes };    // tracking nodes contained in some already found subgraph
    for (int i = 0; i < data->num_subopt_iter; ++i) {
        start_solution = nullptr;                                     // reset to NULL ... Where does the previous one go? Memory leak?
        run_suboptimal_init(&start_solution, &nodes_so_far);          // adapt model for finding next suboptimal, etc.
        solve_successful = model->solve(start_solution);              // try to find next suboptimal solution
        if (!run_suboptimal_windup(solve_successful, &subgraphs, &nodes_so_far, i))      // once no suboptimal subgraph could be found, do not try again ;)
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
            if (getNodeById(data->original_graph, data->nodeid, nodeid, &node))
                if (data->receptors->find(node) != data->receptors->end())
                    subgraph.receptors.insert(nodeid);
        }
    if (data->terminals)
        for (auto nodeid : subgraph.nodes) {
            Node node;
            if (getNodeById(data->original_graph, data->nodeid, nodeid, &node))
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

// Every algorithm type (ie absolute vs average) need its own start heuristic
// TODO: After implementing drgnt in terms of grbfrc::FMILP just branch on Data parameter
// TODO: even more ideally, get rid of specialization altogether? ...

template <> inline
void DeregnetFinder<GRBModel, DrgntData>::find_start_solution(std::pair<Node,
                                                              std::set<Node>>** start_solution) {
    std::function<bool(double, double)> cmp;
    // Adjust comparison function based on whether to minimize or maximize
    if (data->model_sense == "min") cmp = std::greater<double>();
    else cmp = std::less<double>();
    // Try to find start solution and when successful register it
    StartHeuristic heuristic(data->graph,
                             data->score,
                             data->root,
                             data->exclude,
                             data->receptors,
                             cmp,
                             data->size);
    if (heuristic.run())
        *start_solution = heuristic.getStartSolution();
}


template <> inline
void DeregnetFinder<FMILP, AvgdrgntData>::find_start_solution(std::pair<Node,
                                                              std::set<Node>>** start_solution) {


    std::function<bool(double, double)> cmp;
    // Adjust comparison function based on whether to minimize or maximize
    if (data->model_sense == "min") cmp = std::greater<double>();
    else cmp = std::less<double>();
    // Try to find start solution and when successful register it
    AvgStartHeuristic heuristic(data->graph,
                                data->score,
                                data->root,
                                data->exclude,
                                data->receptors,
                                cmp,
                                data->min_size,
                                data->max_size);
    if (heuristic.run())
        *start_solution = heuristic.getStartSolution();
}


// suboptimal heuristics need more information (concerning already found optimal and suboptimal subgraphs) than
// the heuristics for finding optimal subgraphs.

template <> inline
void DeregnetFinder<GRBModel, DrgntData>::find_suboptimal_start_solution(std::pair<Node, std::set<Node>>** start_solution,
                                                                         std::set<std::string>* nodes_so_far) {
    std::function<bool(double, double)> cmp;
    // Adjust comparison function based on whether to minimize or maximize
    if (data->model_sense == "min") cmp = std::greater<double>();
    else cmp = std::less<double>();
    // Try to find start solution and when successful register it
    SuboptimalStartHeuristic heuristic(data->graph,
                                       data->score,
                                       data->root,
                                       data->exclude,
                                       data->receptors,
                                       cmp,
                                       data->nodeid,
                                       nodes_so_far,
                                       data->max_overlap,
                                       data->size);
    if (heuristic.run())
        *start_solution = heuristic.getStartSolution();
}


template <> inline
void DeregnetFinder<FMILP, AvgdrgntData>::find_suboptimal_start_solution(std::pair<Node, std::set<Node>>** start_solution,
                                                                         std::set<std::string>* nodes_so_far) {

    std::function<bool(double, double)> cmp;
    // Adjust comparison function based on whether to minimize or maximize
    if (data->model_sense == "min") cmp = std::greater<double>();
    else cmp = std::less<double>();
    // Try to find start solution and when successful register it
    AvgSuboptimalStartHeuristic heuristic(data->graph,
                                          data->score,
                                          data->root,
                                          data->exclude,
                                          data->receptors,
                                          cmp,
                                          data->nodeid,
                                          nodes_so_far,
                                          data->max_overlap,
                                          data->min_size,
                                          data->max_size);
    if (heuristic.run())
        *start_solution = heuristic.getStartSolution();

}


}    //    namespace deregnet

#endif    //    DEREGNET_FINDER_H
