#ifndef DEREGNET_UTILS_H
#define DEREGNET_UTILS_H

#include <string>
#include <set>
#include <map>

#include <deregnet/usinglemon.h>

namespace deregnet {

// split_string #############################################################

std::set<std::string> split_string(const std::string& str,
                                   const std::string& sep);

// set_union ################################################################

template <typename T>
void set_union(std::set<T>** A, std::set<T> B) {
    for (auto elmt : B) (*A)->insert(elmt);
}

// register_node_set ########################################################

void register_node_set(std::set<std::string>** node_set,
                       std::string comma_separated_list_of_nodes);

// register_node_set_file ###################################################

void register_node_set_file(std::set<std::string>** node_set,
                            std::string path2file);

// ##########################################################################

struct Sbgrph {

    std::string signature;
    std::set<std::pair<std::string, std::string>> edgelist;
    std::set<std::string> nodes;
    std::set<std::string> receptors;
    std::set<std::string> terminals;
    std::string rootid;
    double total_score;
    double avg_score;

    void writeToFile(std::string filepath, bool plain = false);

};

// ##########################################################################

bool getNodeById(Graph* graph, NodeMap<std::string>* nodeid, std::string id, Node* node);

// ##########################################################################

}  // namespace deregnet


#endif    //    DEREGNET_UTILS_H
