#include <iostream>
#include <map>
#include <set>
#include <utility>

#include <gurobi_c++.h>

#include <lemon/concepts/digraph.h>
#include <lemon/smart_graph.h>
#include <lemon/connectivity.h>
#include <lemon/adaptors.h>

#include <deregnet/usinglemon.h>
#include <deregnet/SbgrphModel.h>

using namespace std;

namespace deregnet {

SbgrphModel::SbgrphModel(Graph* graphp,
                         NodeMap<double>* scorep,
                         NodeMap<string>* nodeidp,
                         Node* rootp)
 : graph { graphp },
   score { scorep },
   nodeid { nodeidp },
   root { rootp }
{ }


void SbgrphModel::createVariables() {
    if (root)
        createVariablesRoot();
    else
        createVariablesNoRoot();
    ilp.update();
}

void SbgrphModel::createVariablesRoot() {
    for (NodeIt v(*graph); v != INVALID; ++v)
        x[v] = ilp.addVar(0.0, 1.0, (*score)[v], GRB_BINARY);
}

void SbgrphModel::createVariablesNoRoot() {
    y = new std::map<Node, GRBVar>;
    for (NodeIt v(*graph); v != INVALID; ++v) {
        x[v] = ilp.addVar(0.0, 1.0, (*score)[v], GRB_BINARY);
        (*y)[v] = ilp.addVar(0.0, 1.0, 0.0, GRB_BINARY);
    }
}

void SbgrphModel::addBaseConstraints(int size) {
    if (root)
        addBaseConstraintsRoot(size);
    else
        addBaseConstraintsNoRoot(size);
    ilp.update();
}

void SbgrphModel::addBaseConstraintsRoot(int size) {
    GRBLinExpr subgraph_size_lhs;
    for (NodeIt v(*graph); v != INVALID; ++v) {
        subgraph_size_lhs += x[v];
        GRBLinExpr parent_sum;
        for (InArcIt a(*graph,v); a != INVALID; ++a)
            parent_sum += x[graph->source(a)];
        ilp.addConstr(x[v] - parent_sum <= 0);
    }
    ilp.addConstr(subgraph_size_lhs == size);
}

void SbgrphModel::addBaseConstraintsNoRoot(int size) {
    GRBLinExpr subgraph_size_lhs;
    GRBLinExpr exactly_one_root_lhs;
    for (NodeIt v(*graph); v != INVALID; ++v) {
        ilp.addConstr((*y)[v] <= x[v]);
        subgraph_size_lhs += x[v];
        exactly_one_root_lhs += (*y)[v];
        GRBLinExpr parent_sum;
        for (InArcIt a(*graph,v); a != INVALID; ++a)
            parent_sum += x[graph->source(a)];
        ilp.addConstr(x[v] - (*y)[v] - parent_sum <= 0);
    }
    ilp.addConstr(subgraph_size_lhs == size);
    ilp.addConstr(exactly_one_root_lhs == 1);
}

void SbgrphModel::addIncludeConstraints(std::set<Node>* include) {
    for (Node v : *include)
        ilp.addConstr(x[v] == 1);
    ilp.update();
}

void SbgrphModel::addExcludeConstraints(std::set<Node>* exclude) {
    for (Node v : *exclude)
        ilp.addConstr(x[v] == 0);
    ilp.update();
}

void SbgrphModel::addReceptorConstraints(std::set<Node>* receptors) {
    if (!root)
        for (NodeIt v(*graph); v != INVALID; ++v)
            if (receptors->find(v) == receptors->end())
                ilp.addConstr((*y)[v] == 0);
    ilp.update();
}

void SbgrphModel::addTerminalConstraints(std::set<Node>* terminals) {

}

bool SbgrphModel::solve(bool start_heuristic,
                        double* time_limit,
                        double* gap_cut,
                        std::string model_sense) {
    if (model_sense == "min")
        ilp.set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);
    else
        ilp.set(GRB_IntAttr_ModelSense, GRB_MAXIMIZE);

    ilp.getEnv().set(GRB_IntParam_LazyConstraints, 1);

    if (root) {
        LazyConstraintCallbackRoot callback(&x, graph, root);
        ilp.setCallback(&callback);
        ilp.update();
        ilp.optimize();
    }
    else {
        LazyConstraintCallbackNoRoot callback(&x, y, graph);
        ilp.setCallback(&callback);
        ilp.update();
        ilp.optimize();
    }

    int status = ilp.get(GRB_IntAttr_Status);
    if (status == GRB_OPTIMAL || status == GRB_INTERRUPTED || status == GRB_TIME_LIMIT)
        return true;
    else
        return false;
}

void SbgrphModel::addSuboptimalityConstraint(std::set<std::string> nodes_so_far,
                                             double max_overlap,
                                             int size) {
    GRBLinExpr nodes_so_far_expr;
    for (NodeIt v(*graph); v != INVALID; ++v) {
        if (nodes_so_far.find((*nodeid)[v]) != nodes_so_far.end())
            nodes_so_far_expr += x[v];
    }
    ilp.addConstr(nodes_so_far_expr <= max_overlap * size);
}

SbgrphModel::Solution SbgrphModel::getCurrentSolution() {
    Solution solution;
    set<Node> nodes;
    string rootid;
    for (NodeIt v(*graph); v != INVALID; ++v) {
        if (x[v].get(GRB_DoubleAttr_X) >= 0.98) {
            nodes.insert(v);
            solution.nodes.insert((*nodeid)[v]);
            if ((*y)[v].get(GRB_DoubleAttr_X) >= 0.98)
                rootid = (*nodeid)[v];
        }
    }
    for (auto v : nodes) {
        std::string source = (*nodeid)[v];
        for (OutArcIt a(*graph, v); a != INVALID; ++a) {
            std::string target { (*nodeid)[graph->target(a)] };
            if (solution.nodes.find(target) != solution.nodes.end())
                solution.edgelist.insert( std::make_pair(source, target));
        }
    }
    solution.rootid = rootid;
    solution.total_score = ilp.get(GRB_DoubleAttr_ObjVal);
    solution.avg_score = solution.total_score / solution.nodes.size();
    return solution;
}

// LazyConstraintCallback

SbgrphModel::LazyConstraintCallback::LazyConstraintCallback(std::map<Node, GRBVar>* xx,
                                                            Graph* xgraph,
                                                            Node* xroot)
 : x { xx },
   graph { xgraph },
   root { xroot }
{ }

void SbgrphModel::LazyConstraintCallback::check_and_set_lazy_constr(const int num_components,
                                                                    const InducedSubgraph::NodeMap<int>& component_map) {
    for(int k = 0; k < num_components; k++) {
        std::set<Node> component;
        std::set<Node> parents;
        std::set<Node> global_parents;
        get_component_nodes(component_map, component, k);
        bool is_not_root_component = !is_root_component(component);
        if (is_not_root_component) {
            get_parents(component, parents, global_parents);
            if (parents.size() == 0)
                set_lazy_constraint(component, global_parents);
        }
    }
}

void SbgrphModel::LazyConstraintCallback::callback() {
    try {
        if (where == GRB_CB_MIPSOL) { // new incumbent found
            NodeFilter solution_filter(*graph);
            selected_nodes = {};
            get_solution_filter(solution_filter);
            InducedSubgraph current_subgraph(*graph, solution_filter);
            InducedSubgraph::NodeMap<int> component_map(current_subgraph);
            int num_components { lemon::stronglyConnectedComponents(current_subgraph, component_map) };
            if (num_components > 1)
                check_and_set_lazy_constr(num_components, component_map);
        }
    }
    catch (GRBException e) {
        std::cout << "Error number: " << e.getErrorCode() << std::endl;
        std::cout << e.getMessage() << std::endl;
    }
    catch (...) {
        std::cout << "Error during callback" << std::endl;
    }
}

void SbgrphModel::LazyConstraintCallback::get_component_nodes(const InducedSubgraph::NodeMap<int>& component_map,
                                                              std::set<Node>& component,
                                                              const int k) {
    for (auto v : selected_nodes)
        if (k == component_map[v])
            component.insert(v);
}

bool SbgrphModel::LazyConstraintCallback::is_root_component(const std::set<Node>& component) {
    bool is_root_component = false;
    for (auto v : component)
        if (v == *root) {
            is_root_component = true;
            break;
        }
    return is_root_component;
}

void SbgrphModel::LazyConstraintCallback::get_parents(const std::set<Node>& component,
                                                      std::set<Node>& parents,
                                                      std::set<Node>& global_parents) {
    for (auto v : component) {
        for (InArcIt a(*graph, v); a != INVALID; ++a) {
            Node u = graph->source(a);
            if (component.find(u) == component.end()) {
                global_parents.insert(u);
                if (selected_nodes.find(u) != selected_nodes.end())
                    parents.insert(u);
            }
        }
    }
}


// LazyConstraintCallbackRoot implementation

void SbgrphModel::LazyConstraintCallbackRoot::get_solution_filter(NodeFilter& solution_filter) {
    for (NodeIt v(*graph); v != INVALID; ++v) {
        double current_x = getSolution((*x)[v]);
        if (current_x >= 0.98) {
            solution_filter[v] = true;
            selected_nodes.insert(v);
        }
        else
            solution_filter[v] = false;
    }
}

void SbgrphModel::LazyConstraintCallbackRoot::set_lazy_constraint(const std::set<Node>& component,
                                                                  const std::set<Node>& global_parents) {
    GRBLinExpr component_expr;
    GRBLinExpr global_parents_expr;
    for (auto v : component)
        component_expr += (*x)[v];
    for (auto v : global_parents)
        global_parents_expr += (*x)[v];
    GRBLinExpr lazy_constr_lhs = component_expr - global_parents_expr;
    int comp_size = component.size();
    addLazy(lazy_constr_lhs <= comp_size - 1);
}


// LazyConstraintCallbackNoRoot implementation

SbgrphModel::LazyConstraintCallbackNoRoot::LazyConstraintCallbackNoRoot(std::map<Node, GRBVar>* xx,
                                                                        std::map<Node, GRBVar>* yy,
                                                                        Graph* xgraph)
 : LazyConstraintCallback { xx, xgraph, nullptr },
   y { yy }
{ }



void SbgrphModel::LazyConstraintCallbackNoRoot::get_solution_filter(NodeFilter& solution_filter) {
    for (NodeIt v(*graph); v != INVALID; ++v) {
        double current_x = getSolution((*x)[v]);
        if (current_x >= 0.98) {
            solution_filter[v] = true;
            selected_nodes.insert(v);
        }
        else
            solution_filter[v] = false;
        double current_y = getSolution((*y)[v]);
        if (current_y >= 0.98)
            root = &v;
    }
}


void SbgrphModel::LazyConstraintCallbackNoRoot::set_lazy_constraint(const std::set<Node>& component,
                                                                    const std::set<Node>& global_parents) {
    GRBLinExpr component_expr;
    GRBLinExpr global_parents_expr;
    for (auto v : component)                      // void build_component_expr(GRBLinExpr& expr); ...
        component_expr += (*x)[v] - (*y)[v];
    for (auto v : global_parents)
        global_parents_expr += (*x)[v];
    GRBLinExpr lazy_constr_lhs = component_expr - global_parents_expr;
    int comp_size = component.size();
    addLazy(lazy_constr_lhs <= comp_size - 1);
}

}   // namespace deregnet
