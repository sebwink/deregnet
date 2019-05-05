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

#ifndef DEREGNET_UTILS_H
#define DEREGNET_UTILS_H

#include <string>
#include <set>
#include <map>

#include <deregnet/usinglemon.h>

namespace deregnet {

// split_string #############################################################

/**
 * @brief Split a string into token based on some seperation token
 *
 * @param str String to split
 * @param sep Seperation token  
 *
 * @return The set of corresponding tokens.
 */
std::set<std::string> split_string(const std::string& str,
                                   const std::string& sep);

// set_union ################################################################

/**
 * @brief Inplace set union
 *
 * @tparam T Type of set elements
 * @param[in,out] A First argument of union, union of A and B
 * @param B[in] Second argument of union
 */
template <typename T>
void set_union(std::set<T>** A, const std::set<T>& B) {
    for (auto elmt : B) (*A)->insert(elmt);
}

// register_node_set ########################################################

/**
 * @brief Parse a comma-seperated list of token into a set containing the tokens
 *
 * @param[out] node_set Resulting token set
 * @param[in] comma_separated_list_of_nodes Comma-seperated string of tokens
 */
void register_node_set(std::set<std::string>** node_set,
                       std::string comma_separated_list_of_nodes);

// register_node_set_file ###################################################

/**
 * @brief Parse contents of token file into a set containing the tokens
 *
 * @param[out] node_set Resulting token set
 * @param[in] path2file Path to file, every line is considered a token
 */
void register_node_set_file(std::set<std::string>** node_set,
                            std::string path2file);

// ##########################################################################

/**
 * @brief Representation of a solution to one of the DeRegNet problems 
 */
struct Solution {
    std::set<std::string> nodes;   ///< Nodes in the solution (node id's)
    std::set<std::pair<std::string, std::string>> edgelist;   ///< edgelist (node id --> node id)
    std::string rootid;   ///< Node id of the root node
    double total_score;   ///< Total score, ie \f$\displaystyle \sum_{v \in V_{subgraph}} s_v \f$
    double avg_score;     ///< Total average score, ie total_score / card(nodes)
};

/**
 * @brief Representation of a subgraph corresponding to a Solution
 */
struct Subgraph {

    std::string signature;
    std::set<std::pair<std::string, std::string>> edgelist;
    std::set<std::string> nodes;
    std::set<std::string> receptors;
    std::set<std::string> terminals;
    std::string rootid;
    double total_score;
    double avg_score;

    /**
     * @brief Write to file in SIF format. 
     *
     * @param filepath Path of file
     * @param plain Whether to use plain SIF format or add a metadata header
     */
    void writeToFile(std::string filepath, bool plain = false);

};

// ##########################################################################

/**
 * @brief 
 *
 * @param graph
 * @param nodeid
 * @param id
 * @param node
 *
 * @return 
 */
bool getNodeById(Graph* graph, NodeMap<std::string>* nodeid, std::string id, Node* node);

// ##########################################################################

/**
 * @brief 
 *
 * @param graph
 * @param node_set
 * @param nodeid
 * @param path2file
 */
void write2sif(Graph* graph, std::set<Node>& node_set, NodeMap<std::string>* nodeid, std::string path2file);

}  // namespace deregnet


#endif    //    DEREGNET_UTILS_H
