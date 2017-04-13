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

#include <cmath>
#include <iostream>
#include <fstream>
#include <set>
#include <string>
#include <map>
#include <utility>
#include <algorithm>
#include <functional>

#include <lemon/lgf_reader.h>

#include <gurobi_c++.h>

#include <deregnet/utils.h>
#include <deregnet/usinglemon.h>
#include <deregnet/DeregnetData.h>

using namespace std;

namespace deregnet {

void DeregnetData::read_graph(string* pathToLgf) {
    if (!pathToLgf) {
        cerr << "No lgf specified." << endl;
        exit(INPUT_ERROR);
    }
    graph = new Graph();
    nodeid = new NodeMap<string>(*graph);
    try {
        digraphReader(*graph, *pathToLgf)
        .nodeMap("id", *nodeid)
        .run();
    }
    catch ( ... ) {
        cerr << "Could not read lgf file " << pathToLgf << "." << endl;
        exit(INPUT_ERROR);
    }
    if (!receptor_as_root)
        reverse_graph();
}

void DeregnetData::read_score(string* pathToTsv, bool take_abs) {
    if (!graph) {
        cerr << "Error: call AvgDeregnetData::read_graph() first." << endl;
        exit(INPUT_ERROR);
    }
    if (!pathToTsv) {
        cerr << "No score file specified." << endl;
        exit(INPUT_ERROR);
    }
    map<string, double> score_map;
    try {
        ifstream tsv;
        tsv.open( *pathToTsv );
        string key;
        double value;
        while ( !tsv.eof() ) {
              tsv >> key >> value;
              std::cout <<  key << " : " << value << std::endl;
              score_map[key] = value;
        }
        score = new NodeMap<double>(*graph);
        for (NodeIt v(*graph); v != INVALID; ++v) {
            (*score)[v] = score_map[(*nodeid)[v]];
            if (take_abs)
                (*score)[v] = std::abs((*score)[v]);
        }
        tsv.close();
    }
    catch ( ... ) {
      cerr << "Could not read score file " << pathToTsv << "." << endl;
      exit(INPUT_ERROR);
    }
}

void DeregnetData::get_node_set(std::set<Node>** node_set,
                                      std::set<std::string>* node_ids) {
    if (!graph) {
        cerr << "Error: call DeregnetData::read_graph() first." << endl;
        exit(INPUT_ERROR);
    }
    if (node_ids) {
        if (!(*node_set))
            *node_set = new set<Node>();
        Node next;
        for (auto id : *node_ids)
            if ( getNodeById(graph, nodeid, id, &next) )
                (*node_set)->insert(next);
    }
}

void DeregnetData::get_root(std::string* root_id) {
    if (root_id) {
        root = new Node();
        if (!getNodeById(graph, nodeid, *root_id, root)) {
            cerr << "Root seems not to be in the graph." << endl;
            exit(INPUT_ERROR);
        }
    }
}

void DeregnetData::reverse_graph() {
    map<Node,Node> node_correspondence;
    revgraph = new Graph();
    revnodeid = new NodeMap<string>(*revgraph);
    for (NodeIt v(*graph); v != INVALID; ++v) {
        Node revv = revgraph->addNode();
        (*revnodeid)[revv] = (*nodeid)[v];
        node_correspondence[v] = revv;
    }
    for (NodeIt v(*graph); v != INVALID; ++v) {
        Node revv = node_correspondence[v];
        for (InArcIt a(*graph, v); a != INVALID; ++a) {
            Node s = graph->source(a);
            Node revs = node_correspondence[s];
            revgraph->addArc(revv,revs);
        }
    }
    original_graph = graph;
    original_nodeid = nodeid;
    graph = revgraph;
    nodeid = revnodeid;
    swap(receptors, terminals);
}

}    //    namespace deregnet
