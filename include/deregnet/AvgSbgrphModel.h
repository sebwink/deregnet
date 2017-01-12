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

#ifndef AVG_SBGRPH_MODEL_H
#define AVG_SBGRPH_MODEL_H

#include <string>
#include <map>
#include <set>
#include <utility>

#include <grbfrc/FMILP.h>
#include <deregnet/usinglemon.h>

namespace deregnet {

class AvgSbgrphModel {

  public:

    typedef struct {
        std::set<std::string> nodes;
        std::set<std::pair<std::string,std::string>> edgelist;
        std::string rootid;
        double total_score;
        double avg_score;
    } Solution;

  private:

    GRBEnv env;
    grbfrc::FMILP filp { grbfrc::FMILP(env) };
    std::map<Node, GRBVar> x;
    Graph* graph;
    NodeMap<double>* score;
    NodeMap<std::string>* nodeid;
    Node* root { nullptr };
    std::map<Node, GRBVar>* y { nullptr };
/*
    class LazyConstraintCallback : public GrbfrcCallback {

    };
*/
  public:

    AvgSbgrphModel(Graph* graphp,
                   NodeMap<double>* scorep,
                   NodeMap<string>* nodeidp,
                   Node* rootp);
    void createVariables();
    void reverseGraph();
    void addBaseConstraints(int min_size,
                            int max_size);
    void addIncludeConstraints(std::set<Node>* include);
    void addExcludeConstraints(std::set<Node>* exclude);
    void addReceptorConstraints(std::set<Node>* receptors);
    void addTerminalConstraints(std::set<Node>* terminals);
    bool solve(grbfrc::Algorithm algorithm,
               bool start_heuristic,
               double* time_limit,
               double* gap_cut,
               std::string model_sense);
    void addSuboptimalityConstraint(std::set<std::string> nodes_so_far,
                                    double max_overlap);
    Solution getCurrentSolution();

  private:

    void createVariablesRoot();
    void createVariablesNoRoot();
    void addBaseConstraintsRoot(int min_size, int max_size);
    void addBaseConstraintsNoRoot(int min_size, int max_size);

};


}   // namespace deregnet

#endif    // AVG_SBGRPH_MODEL_H
