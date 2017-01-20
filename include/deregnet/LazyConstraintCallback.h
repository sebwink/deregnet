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

#ifndef LAZY_CONSTRAINT_CALLBACK_H
#define LAZY_CONSTRAINT_CALLBACK_H

#include <string>
#include <map>
#include <vector>
#include <set>

#include <deregnet/usinglemon.h>

#include <gurobi_c++.h>

namespace deregnet {

class LazyConstraintCallback : public GRBCallback {

  public:

    LazyConstraintCallback(std::map<Node, GRBVar>* xx,
                           Graph* xgraph,
                           Node* xroot,
                           double* xgap_cut);
  protected:

    std::map<Node, GRBVar>* x;
    Graph* graph;
    Node* original_root;
    Node* root;
    std::vector<Node> selected_nodes;
    double* gap_cut;

  protected:

    void callback();

 private:

    virtual void get_solution_filter(NodeFilter& solution_filter) = 0;
    virtual void set_lazy_constraint(const std::set<Node>& component,
                                     const std::set<Node>& global_parents) = 0;
    void check_and_set_lazy_constr(const int num_components,
                                   const InducedSubgraph::NodeMap<int>& component_map);
    void get_component_nodes(const InducedSubgraph::NodeMap<int>& component_map,
                             std::set<Node>& component,
                             const int k);
    bool is_root_component(const std::set<Node>& component);
    void get_parents(const std::set<Node>& component,
                     std::set<Node>& parents,
                     std::set<Node>& global_parents);
};

class LazyConstraintCallbackRoot : public LazyConstraintCallback {

  public:

    using LazyConstraintCallback::LazyConstraintCallback;

  private:

    virtual void get_solution_filter(NodeFilter& solution_filter) override;
    virtual void set_lazy_constraint(const std::set<Node>& component,
                                     const std::set<Node>& global_parents) override;

};

class LazyConstraintCallbackNoRoot : public LazyConstraintCallback {

  public:

    LazyConstraintCallbackNoRoot(std::map<Node, GRBVar>* xx,
                                 std::map<Node, GRBVar>* yy,
                                 Graph* xgraph,
                                 double* xgap_cut);

  private:

    std::map<Node, GRBVar>* y;

    virtual void get_solution_filter(NodeFilter& solution_filter) override;
    virtual void set_lazy_constraint(const std::set<Node>& component,
                                     const std::set<Node>& global_parents) override;

};

}    //    namespace deregnet

#endif    //    LAZY_CONSTRAINT_CALLBACK_H
