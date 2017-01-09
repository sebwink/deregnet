#ifndef SBGRPH_FINDER_H
#define SBGRPH_FINDER_H

#define INPUT_ERROR 1
#define NO_OPTIMAL_SUBGRAPH_FOUND 2

#include <set>
#include <string>
#include <utility>

#include <deregnet/usinglemon.h>
#include <deregnet/utils.h>
#include <deregnet/SbgrphModel.h>

namespace deregnet {


class SbgrphFinder {

  public:

    struct Data {
        Graph* graph { nullptr };
        NodeMap<double>* score { nullptr };
        NodeMap<std::string>* nodeid { nullptr };
        Node* root { nullptr };
        int size { 20 };
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
        Graph* revgraph;
        NodeMap<double>* revscore { nullptr };
        NodeMap<std::string>* revnodeid { nullptr };
        Node* revroot { nullptr };

        void get_reversed_graph();
        void read_graph(std::string* pathToLgf);
        void read_score(std::string* pathToTsv);
        void get_node_set(std::set<Node>* node_set, std::set<std::string>* node_ids);
        void get_root(std::string* root_id);

    };

  private:

    Data data;

    Sbgrph toSbgrph(SbgrphModel::Solution solution, std::string signature);

  public:

    SbgrphFinder(Data data);
    std::vector<Sbgrph> run(bool start_heuristic, std::string model_sense);


};


}   // namespace deregnet

#endif   // SBGRPH_FINDER_H
