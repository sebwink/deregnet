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

void Subgraph::writeToFile(string filepath, bool plain) {
    ofstream sif_file;
    sif_file.open(filepath + "/" + signature + ".sif");
    if (!plain) {
        sif_file << "# root : " << rootid << "\n"
                 << "# total_score : " << total_score << "\n"
                 << "# avg_score : " << avg_score << "\n"
                 << "# receptors : ";
        for (auto receptor: receptors)
            sif_file << receptor << " ";
        sif_file << "\n# terminals : ";
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

void write2sif(Graph* graph, std::set<Node>& node_set, NodeMap<std::string>* nodeid, std::string path2file) {
    ofstream sif_file;
    sif_file.open(path2file);
    for (auto v : node_set)
        for (OutArcIt a(*graph, v); a != INVALID; ++a)
            if (node_set.find(graph->target(a)) != node_set.end())
                sif_file << (*nodeid)[v] << "\tpp\t" << (*nodeid)[graph->target(a)] << "\n";
    sif_file.close();
}

}   //  namespace deregnet
