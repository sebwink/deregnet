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

#ifndef DEREGNET_DATA_H
#define DEREGNET_DATA_H

#define INPUT_ERROR 1
#define NO_OPTIMAL_SUBGRAPH_FOUND 2

#include <set>
#include <string>
#include <utility>

#include <deregnet/utils.h>
#include <deregnet/usinglemon.h>
#include <deregnet/DeregnetModel.h>

namespace deregnet {

class DeregnetData {

      public:

        Graph* graph { nullptr };
        NodeMap<double>* score { nullptr };
        NodeMap<std::string>* nodeid { nullptr };
        Node* root { nullptr };
        std::set<Node>* terminals { nullptr };
        std::set<Node>* receptors { nullptr };
        std::set<Node>* include { nullptr };
        std::set<Node>* exclude { nullptr };
        int num_subopt_iter { 0 };                         /*< number of iterations to find suboptimal subgraphs */
        double max_overlap { 0.0 };                        /*< maximal percentage of overlap to previous subgraphs
                                                               when searching for suboptimal subgraphs             */
        double* time_limit { nullptr };                    /*< time limit of a single solve */
        double* gap_cut { nullptr };                       /*< gap tolerance */

        bool receptor_as_root { true };

     private:

        Graph* revgraph { nullptr };
        NodeMap<std::string>* revnodeid { nullptr };

        Graph* original_graph { nullptr };
        NodeMap<std::string>* original_nodeid { nullptr };

      public:

        void read_graph(std::string* pathToLgf);
        void read_score(std::string* pathToTsv);
        void get_node_set(std::set<Node>** node_set, std::set<std::string>* node_ids);
        void get_root(std::string* root_id);

      private:

        void reverse_graph();

    };

}   // namespace deregnet

#endif   // DEREGNET_DATA_H
