#ifndef AVG_SBGRPH_FINDER_H
#define AVG_SBGRPH_FINDER_H

#define INPUT_ERROR 1
#define NO_OPTIMAL_SUBGRAPH_FOUND 2

#include <set>
#include <string>

#include <deregnet/usinglemon.h>
#include <deregnet/utils.h>
#include <deregnet/AvgSbgrphModel.h>

namespace deregnet {


class AvgSbgrphFinder {

  public:

    typedef struct {
        Graph* graph { nullptr };
        NodeMap<double>* score { nullptr };
        NodeMap<std::string>* nodeid { nullptr };
        Node* root { nullptr };
        int max_size { 50 };                               /*< maximal size of subgraph */
        int min_size { 20 };                               /*< minimal size of subgraph */
        std::set<Node>* terminals { nullptr };
        std::set<Node>* receptors { nullptr };
        std::set<Node>* include { nullptr };
        std::set<Node>* exclude { nullptr };
        int num_subopt_iter { 0 };                         /*< number of iterations to find suboptimal subgraphs */
        double max_overlap { 0.0 };                        /*< maximal percentage of overlap to previous subgraphs
                                                               when searching for suboptimal subgraphs             */
        double* time_limit { nullptr };                    /*< time limit of a single solve */
        double* gap_cut { nullptr };                       /*< gap tolerance, not suitable for algorithm dta */

        void read_graph(std::string* pathToLgf);
        void read_score(std::string* pathToTsv);
        void get_node_set(std::set<Node>* node_set, std::set<std::string>* node_ids);
        void get_root(std::string* root_id);

    } Data;

  private:

    Data data;
    bool receptor_as_root;

    Sbgrph toSbgrph(AvgSbgrphModel::Solution solution, std::string signature);

  public:

    AvgSbgrphFinder(Data data, bool receptor_as_root);
    std::vector<Sbgrph> run(grbfrc::Algorithm algorithm, bool start_heuristic, std::string model_sense);


};


}   // namespace deregnet

#endif   // AVG_SBGRPH_FINDER_H
