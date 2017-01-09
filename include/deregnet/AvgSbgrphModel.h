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
