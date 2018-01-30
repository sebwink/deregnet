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

#ifndef DEREGNET_START_HEURISTIC_H
#define DEREGNET_START_HEURISTIC_H

#include <set>
#include <utility>
#include <functional>

#include <deregnet/usinglemon.h>

namespace deregnet {

/**
 * @brief Abstract base class for greedy start heuristics
 */
class DeregnetStartHeuristic {

  protected:

    Graph* graph;       ///< Underlying graph (in which we try to find subgraphs in)
    NodeMap<double>* score;     ///< Node scores
    Node* root;         ///< Root node (if available)
    Node* original_root;        ///< Original root node (passed on construction)
    std::set<Node>* exclude;    ///< Nodes to exclude to from any subgraph
    std::set<Node>* receptors;  ///< Nodes considered to act as 'receptors'

    std::function<bool(double, double)> cmp;    ///< Comparison function to compare solutions

    std::set<Node>* start_solution { nullptr };     ///< heuristic start solution as a set of nodes (if one is found)
    bool found_node { false };  ///< Internal variable indicating progress of the heuristic
    bool success { false };     ///< Indicator whether a heuristic solution could be found

  public:

    /**
     * @brief Constructor
     *
     * @param xgraph Underlying network (in which we try to find subgraphs in)
     * @param xscore Node scores
     * @param root Root node (if available)
     * @param exclude Nodes to exclude from any subgraph (if available)
     * @param receptors Nodes considered to be 'receptor' nodes
     * @param xcmp Comparison function (to compare solutuions, eg std:less)
     */
    DeregnetStartHeuristic(Graph* xgraph,
                           NodeMap<double>* xscore,
                           Node* root,
                           std::set<Node>* exclude,
                           std::set<Node>* receptors,
                           std::function<bool(double, double)> xcmp);
    
    /**
     * @brief Run the heuristic (ie try to find a heuristic start solution)
     *
     * @return Whether a heuristic solution was found
     */
    bool run();

    /**
     * @brief Retrieve the start solution (if found)
     *
     * @return Tuple containing root node of solution and the solution itself
     */
    std::pair<Node, std::set<Node>>* getStartSolution();

  protected:

    /**
     * @brief Find next node to add to solution
     *
     * @return A node (nullptr if none could be found)
     */
    Node* get_next_node();
    
    /**
     * @brief Update status of 'best' node
     *
     * @param current_best_node Variable capturing the current best node while execution of the algorithm
     * @param node Node to become the new best node
     * @param best Score value of the best node
     */
    void update_best_node(Node** current_best_node, Node* node, double* best);

    /**
     * @brief Find best node to add to solution in every iteration
     *
     * @return The best node to add
     */
    virtual Node* get_best_root() = 0;
    
    /**
     * @brief Decide whether to continue searching further (for other solution than the one already found)
     *
     * @return Whether to search further
     */
    virtual bool search_further() = 0;
    
    /**
     * @brief Decide whether adding a given node to the solution is feasible wrt the given constraints
     *
     * @param node A node in the graph
     *
     * @return Whether adding the node is feasible
     */
    virtual bool feasible_node(Node* node) = 0;

};

}    //    namespace deregnet

#endif    //    DEREGNET_START_HEURISTIC
