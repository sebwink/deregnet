#include <fstream>
#include <string>
#include <set>
#include <map>

#include <deregnet/usinglemon.h>
#include <deregnet/utils.h>

using namespace std;

namespace deregnet {

// split_string #############################################################

set<string> split_string(const string& str,
                         const string& sep) {
    set<string> split;
    string::size_type start { 0 };
    string::size_type last_sep { 0 };
    while (last_sep != string::npos) {
        last_sep = str.find_first_of(sep, start);
        split.insert( str.substr(start, last_sep - start) );
        start = last_sep + 1;
    }
    return split;
}

// register_node_set ########################################################

void register_node_set(set<string>** node_set,
                       string comma_separated_list_of_nodes) {
    if (!(*node_set))
        *node_set = new set<string>( split_string(comma_separated_list_of_nodes, ",") );
    else
        set_union<string>(node_set, split_string(comma_separated_list_of_nodes, ",") );
}

// register_node_set_file ###################################################

void register_node_set_file(set<string>** node_set,
                            string path2file) {
    set<string> nodes_in_file;
    try {
        ifstream file;
        file.open( path2file );
        string nodeid;
        while ( !file.eof() ) {
            file >> nodeid;
            nodes_in_file.insert(nodeid);
        }
        file.close();
    }
    catch ( ... ) {
        cerr << "Could not read file " << path2file << endl;
    }
    if (!(*node_set))
        *node_set = new set<string>(nodes_in_file);
    else
        set_union<string>(node_set, nodes_in_file);
}

// Sbgrph::writeToFile #######################################################

void Sbgrph::writeToFile(string filepath, bool plain) {
    ofstream sif_file;
    sif_file.open(filepath + "/" + signature + ".sif");
    if (!plain) {
        sif_file << "# root : " << rootid << "\n"
                 << "# total_score : " << total_score << "\n"
                 << "# avg_score : " << avg_score << "\n"
                 << "# receptors : ";
        for (auto receptor: receptors)
            sif_file << receptor << " ";
        sif_file << "\n#terminals : ";
        for (auto terminal: terminals)
            sif_file << terminal << " ";
        sif_file << "\n\n";
    }
    for (auto edge : edgelist) {
        sif_file << edge.first
                 << "\tpp\t"
                 << edge.second
                 << "\n";
    }
    sif_file.close();
}

// getNodeById ##############################################################

bool getNodeById(Graph* graph, NodeMap<std::string>* nodeid, std::string id, Node* node) {
    for (NodeIt v(*graph); v != INVALID; ++v)
        if ((*nodeid)[v] == id) {
            *node = v;
            return true;
        }
    return false;
}

// ##########################################################################


}   //  namespace deregnet
