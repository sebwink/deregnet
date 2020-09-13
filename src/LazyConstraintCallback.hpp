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

#ifndef LAZY_CONSTRAINT_CALLBACK_H
#define LAZY_CONSTRAINT_CALLBACK_H

#include <string>
#include <map>
#include <vector>
#include <set>

#include "usinglemon.hpp"

#include <gurobi_c++.h>

namespace deregnet {

/**
 * @brief Base lazy constraint callback class.
 *
 * Implements most of the basic logic concerning lazy constraints.
 * Specializations arise from whether there is a fixed root or not
 * (see LazyConstraintCallbackRoot and LazyConstraintCallbackNoRoot)
 * and/or whether variable groups have to be mapped with respect to
 * reformulations of the original problem (see GrbfrcLazyConstraintCallbackRoot
 * and GrbfrcLazyConstraintCallbackNoRoot).
 *
 */
class LazyConstraintCallback : public GRBCallback {

  public:

    /**
     * @brief Constructor passing data relevant for the callback. 
     *
     * @param xx Mapping from graph nodes to subgraph containment indicator variables
     * @param xgraph Graph one is searching for subgraphs in
     * @param xroot The root node (or null if there is no fixed one)
     * @param xgap_cut The gap cut parameter
     */
    LazyConstraintCallback(std::map<Node, GRBVar>* xx,
                           Graph* xgraph,
                           Node* xroot,
                           double* xgap_cut);
  protected:

    std::map<Node, GRBVar>* x;      ///< Mapping from graph nodes to subgraph containment indicator variables
    Graph* graph;                   ///< The underlying graph we are searching for subgraphs in
    Node* original_root;            ///< Root at beginning of the callback
    Node* root;                     ///< Current root node
    std::vector<Node> selected_nodes;    ///< Current nodes in the subgraph
    double* gap_cut;                ///< Gap cut value

  protected:

    /**
     * @brief The actual callback function passed to the Gurobi optimizer
     *
     * Describe overall callback algorithm
     *
     */
    void callback();

 private:

    /**
     * @brief Figure out the current solution
     *
     * This depends on whether there is a fixed root. Correspondingly, 
     * conrete callback classed must implement this method in a suitable
     * fashion.
     *
     * @param[out] solution_filter Boolean node map indicating nodes in current solution/subgraph
     */
    virtual void get_solution_filter(NodeFilter& solution_filter) = 0;
    
    /**
     * @brief Add lazy constraints to the model.
     *
     * Given a solution \f$G'=(V',E')\f$ and a strong component \f$S \subset V'\f$ such that
     * for all nodes \f$u \in V' \backslash S \f$, if \f$(u,v) \in E'\f$ then \f$v \notin S\f$, we add
     * some lazy constraints to the model. These constraints have a different mathematical
     * form depending on whether the root is fixed or not and hence concrete callback classes
     * must implement this method correspondingly.
     *
     * @param component An strong component of the solution not reachable from other nodes of the solution
     * @param global_parents All nodes not in the component which have a an edge towards a node from the component
     */
    virtual void set_lazy_constraint(const std::set<Node>& component,
                                     const std::set<Node>& global_parents) = 0;
    
    /**
     * @brief Check whether strong components violate any constraints and if so, add them
     *
     * @param num_components Number of strong components in the current subgraph
     * @param component_map Data structure encoding connectivity of the subgraph, mapping nodes to component id's
     */
    void check_and_set_lazy_constr(const int num_components,
                                   const InducedSubgraph::NodeMap<int>& component_map);
    
    /**
     * @brief Gets the nodes belonging to a strong component 
     *
     * @param[in] component_map Data structure encoding connectivity of the subgraph, mapping nodes to components id's
     * @param[out] component Nodes in the component
     * @param[in] k Component id
     *
     * @return Whether the component consists of more than one node 
     */
    bool get_component_nodes(const InducedSubgraph::NodeMap<int>& component_map,
                             std::set<Node>& component,
                             const int k);
    
    /**
     * @brief Check whether a component contains the root (no constraint violation possible in that case)
     *
     * @param component A strong component of the current subgraph
     *
     * @return Whether the component contains the root
     */
    bool is_root_component(const std::set<Node>& component);
    
    /**
     * @brief Find in-neighbors of a strong component within the subgraph and globally 
     *
     * @param[in] component A strong component of the current subgraph
     * @param[out] parents Subgraph nodes (not in component) having edges connecting to a component node
     * @param[out] global_parents Any nodes (not in component) having edges connecting to a component node
     */
    void get_parents(const std::set<Node>& component,
                     std::set<Node>& parents,
                     std::set<Node>& global_parents);
};

/**
 * @brief Implementation of the callback when there is a fixed root
 */
class LazyConstraintCallbackRoot : public LazyConstraintCallback {

  public:

    /**
     * @brief Inherit constructor from LazyConstraintCallback
     */
    using LazyConstraintCallback::LazyConstraintCallback;

  private:

    /**
     * @brief Figure out the nodes in the current subgraph
     *
     * @param[out] solution_filter Node mapping indicating whether a node belongs to the current solution
     */
    virtual void get_solution_filter(NodeFilter& solution_filter) override;
    
    /**
     * @brief Add lazy constraints to model.
     *
     * Given a solution \f$G'=(V',E') \subset G=(V,E)\f$ and a strong component \f$S \subset V'\f$ such that
     * for all nodes \f$u \in V' \backslash S \f$, if \f$(u,v) \in E'\f$ then \f$v \notin S\f$, we add
     * the following lazy constraint to the model:
     * \f[
     *      \begin{equation}
     *          \displaystyle \sum_{v \in S} x_v - \displaystyle \sum_{u \in N^{-}(S)} x_u \leq |S| - 1
     *      \end{equation}
     * \f]
     * Here, \f$N^{-}(S) = \{v \in V \backslash S: \exists u \in S: (v,u) \in E \} \f$ are the global (incoming)
     * parents of the component \f$S\f$.
     *
     * @param component An strong component of the solution not reachable from other nodes of the solution
     * @param global_parents All nodes not in the component which have a an edge towards a node from the component
     */
    virtual void set_lazy_constraint(const std::set<Node>& component,
                                     const std::set<Node>& global_parents) override;

};

/**
 * @brief Implementation of the callback when there is no fixed root
 */
class LazyConstraintCallbackNoRoot : public LazyConstraintCallback {

  public:

    /**
     * @brief 
     *
     * @param xx
     * @param yy
     * @param xgraph
     * @param xgap_cut
     */
    LazyConstraintCallbackNoRoot(std::map<Node, GRBVar>* xx,
                                 std::map<Node, GRBVar>* yy,
                                 Graph* xgraph,
                                 double* xgap_cut);

  private:

    std::map<Node, GRBVar>* y;      ///< Mapping from graph nodes to root indicator variables

    /**
     * @brief Figure out the current nodes in the subgraph and the identity of the root 
     *
     * @param[out] solution_filter Node mapping indicating whether a node belongs to the current solution
     */
    virtual void get_solution_filter(NodeFilter& solution_filter) override;
    
    /**
     * @brief Add lazy constraints to model.
     *
     * Given a solution \f$G'=(V',E') \subset G=(V,E)\f$ and a strong component \f$S \subset V'\f$ such that
     * for all nodes \f$u \in V' \backslash S \f$, if \f$(u,v) \in E'\f$ then \f$v \notin S\f$, we add
     * the following lazy constraints to the model:
     * \f[
     *      \begin{equation}
     *          \displaystyle \sum_{v \in S} (x_v - y_v) - \displaystyle \sum_{u \in N^{-}(S)} x_u \leq |S| - 1
     *      \end{equation}
     * \f]
     * Here, \f$N^{-}(S) = \{v \in V \backslash S: \exists u \in S: (v,u) \in E \} \f$ are the global (incoming)
     * parents of the component \f$S\f$.
     * 
     * @param component An strong component of the solution not reachable from other nodes of the solution
     * @param global_parents All nodes not in the component which have a an edge towards a node from the component
     */
    virtual void set_lazy_constraint(const std::set<Node>& component,
                                     const std::set<Node>& global_parents) override;

};

}    //    namespace deregnet

#endif    //    LAZY_CONSTRAINT_CALLBACK_H
