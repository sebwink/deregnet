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

#ifndef DEREGNET_UTILS_H
#define DEREGNET_UTILS_H

#include <string>
#include <set>
#include <map>

#include <deregnet/usinglemon.h>

namespace deregnet {

// split_string #############################################################

std::set<std::string> split_string(const std::string& str,
                                   const std::string& sep);

// set_union ################################################################

template <typename T>
void set_union(std::set<T>** A, std::set<T> B) {
    for (auto elmt : B) (*A)->insert(elmt);
}

// register_node_set ########################################################

void register_node_set(std::set<std::string>** node_set,
                       std::string comma_separated_list_of_nodes);

// register_node_set_file ###################################################

void register_node_set_file(std::set<std::string>** node_set,
                            std::string path2file);

// ##########################################################################

struct Solution {
    std::set<std::string> nodes;
    std::set<std::pair<std::string, std::string>> edgelist;
    std::string rootid;
    double total_score;
    double avg_score;
};

struct Subgraph {

    std::string signature;
    std::set<std::pair<std::string, std::string>> edgelist;
    std::set<std::string> nodes;
    std::set<std::string> receptors;
    std::set<std::string> terminals;
    std::string rootid;
    double total_score;
    double avg_score;

    void writeToFile(std::string filepath, bool plain = false);

};

// ##########################################################################

bool getNodeById(Graph* graph, NodeMap<std::string>* nodeid, std::string id, Node* node);

// ##########################################################################

void write2sif(Graph* graph, std::set<Node>& node_set, NodeMap<std::string>* nodeid, std::string path2file);

}  // namespace deregnet


#endif    //    DEREGNET_UTILS_H
