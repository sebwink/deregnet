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

#ifndef AVG_SUBOPTIMAL_START_HEURISTIC_H
#define AVG_SUBOPTIMAL_START_HEURISTIC_H

#include <string>
#include <set>
#include <utility>

#include <deregnet/usinglemon.h>
#include <deregnet/DeregnetStartHeuristic.h>

namespace deregnet {

class AvgSuboptimalStartHeuristic : public DeregnetStartHeuristic {

  private:

    NodeMap<std::string>* nodeid;
    std::set<std::string>* nodes_so_far;
    double max_overlap;
    int min_size;
    int max_size;
    int current_overlapping_nodes { 0 };
    double current_avg { 0.0 };

  public:

    AvgSuboptimalStartHeuristic(Graph* xgraph,
                                NodeMap<double>* xscore,
                                Node* root,
                                std::set<Node>* exclude,
                                std::set<Node>* receptors,
                                std::function<bool(double, double)> xcmp,
                                NodeMap<std::string>* xnodeid,
                                std::set<std::string>* xnodes_sof_far,
                                double xmax_overlap,
                                int xmin_size,
                                int xmax_size);

  private:

    bool is_overlap_node(Node* node);
    virtual Node* get_best_root() override;
    virtual bool search_further() override;
    virtual bool feasible_node(Node* node) override;

    double update_avg(Node* node);

};

}    //    namespace deregnet

#endif    //    AVG_SUBOPTIMAL_START_HEURISTIC
