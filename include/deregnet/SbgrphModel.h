#ifndef SBGRPH_MODEL_H
#define SBGRPH_MODEL_H

#include <string>
#include <map>
#include <set>
#include <utility>

#include <gurobi_c++.h>

#include <deregnet/usinglemon.h>
#include <lemon/concepts/digraph.h>

namespace deregnet {

class SbgrphModel {

  public:

    struct Solution {
        std::set<std::string> nodes;
        std::set<std::pair<std::string, std::string>> edgelist;
        std::string rootid;
        double total_score;
        double avg_score;
    };

  private:

    GRBEnv env;
    GRBModel ilp { GRBModel(env) };
    std::map<Node, GRBVar> x;
    Graph* graph;
    NodeMap<double>* score;
    NodeMap<std::string>* nodeid;
    Node* root { nullptr };
    std::map<Node, GRBVar>* y { nullptr };

    class LazyConstraintCallback : public GRBCallback {

      public:

        LazyConstraintCallback(std::map<Node, GRBVar>* xx,
                               Graph* xgraph,
                               Node* xroot);
      protected:

        std::map<Node, GRBVar>* x;
        Graph* graph;
        Node* root;
        std::set<Node> selected_nodes;

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
                                     Graph* xgraph);

      private:

        std::map<Node, GRBVar>* y;

        virtual void get_solution_filter(NodeFilter& solution_filter) override;
        virtual void set_lazy_constraint(const std::set<Node>& component,
                                         const std::set<Node>& global_parents) override;

    };

  public:

    SbgrphModel(Graph* graphp,
                NodeMap<double>* scorep,
                NodeMap<string>* nodeidp,
                Node* rootp);
    void createVariables();
    void addBaseConstraints(int size);
    void addIncludeConstraints(std::set<Node>* include);
    void addExcludeConstraints(std::set<Node>* exclude);
    void addReceptorConstraints(std::set<Node>* receptors);
    void addTerminalConstraints(std::set<Node>* terminals);
    bool solve(bool start_heuristic,
               double* time_limit,
               double* gap_cut,
               std::string model_sense);
    void addSuboptimalityConstraint(std::set<std::string> nodes_so_far,
                                    double max_overlap,
                                    int size);
    Solution getCurrentSolution();

  private:

    void createVariablesRoot();
    void createVariablesNoRoot();
    void addBaseConstraintsRoot(int size);
    void addBaseConstraintsNoRoot(int size);

};


}   // namespace deregnet

#endif    // SBGRPH_MODEL_H
