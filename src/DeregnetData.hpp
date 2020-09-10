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

#ifndef DEREGNET_DATA_H
#define DEREGNET_DATA_H

#define INPUT_ERROR 1
#define NO_OPTIMAL_SUBGRAPH_FOUND 2

#include <utility>

#include "utils.hpp"

namespace deregnet {

/**
 * @brief Represents all shared parameters and data needed to run a DeRegNet algorithm.
 *
 * Captures all data shared among all DeRegNet algorithms and allows
 * to read relevant data from files and carry out certain pre-processing
 * steps, like reversing the input graph's edge orientations.
 *
 */
class DeregnetData {

  public:

    Graph* graph { nullptr };                          ///< Actual graph used to run the algorithms
    Graph* original_graph { nullptr };                 ///< Original input graph
    NodeMap<double>* score { nullptr };                ///< Node scores 
    NodeMap<std::string>* nodeid { nullptr };          ///< Node id's
    Node* root { nullptr };                            ///< Root node (null if the algorithm is supposed to find the root)
    std::set<Node>* terminals { nullptr };             ///< Nodes classified as 'terminals'
    std::set<Node>* receptors { nullptr };             ///< Nodes classified as 'receptors'
    std::set<Node>* include { nullptr };               ///< Nodes included in any subgraph a priori
    std::set<Node>* exclude { nullptr };               ///< Nodes excluded from any subgraph a priori
    int num_subopt_iter { 0 };                         ///< Number of iterations to find suboptimal subgraphs 
    double max_overlap { 0.0 };                        ///< Maximal percentage of overlap to previous subgraphs when searching for suboptimal subgraphs
    double* time_limit { nullptr };                    ///< Time limit of a single solve
    double* gap_cut { nullptr };                       ///< Gap tolerance 
    bool receptor_as_root { true };                    ///< Whether to orient the subgraphs such that the root acts as 'receptor'
    std::string model_sense { "max" };                 ///< Whether to maximize or minimize
    bool start_heuristic { true };                     ///< Whether to run the greedy start heuristic

  private:

    NodeMap<std::string>* revnodeid { nullptr };
    Graph* revgraph { nullptr };

    NodeMap<std::string>* original_nodeid { nullptr };

  public:

    /**
     * @brief Read graph from LGF (Lemon Graph Format) file
     *
     * http://lemon.cs.elte.hu/pub/tutorial/a00018.html
     *
     * @param pathToLgf Path to LGF file
     */
    void read_graph(std::string* pathToLgf);
    
    /**
     * @brief Read scores from file. 
     *
     * @param pathToTsv Path to file with scores: two columns tab-seperated, first column: node id's
     * @param take_abs Whether to take the absolute values of the parsed scores as actual score
     */
    void read_score(std::string* pathToTsv, bool take_abs);
   
    /**
     * @brief Translate a set of node id's into a set of (Lemon) nodes 
     *
     * @param[out] node_set Resulting node set
     * @param[in] node_ids Set of node id's
     */
    void get_node_set(std::set<Node>** node_set, std::set<std::string>* node_ids);
    
    /**
     * @brief Get (Lemon) node corresponding to the id of the root 
     *
     * Updates DeregnetData.root
     *
     * @param root_id Node id of the root node (null if none specified)
     */
    void get_root(std::string* root_id);

  private:

    /**
     * @brief Reverse orientation of edges in original input graph
     */
    void reverse_graph();

};

}   // namespace deregnet

#endif   // DEREGNET_DATA_H
