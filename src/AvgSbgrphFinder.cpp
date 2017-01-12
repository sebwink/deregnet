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

#include <iostream>
#include <fstream>
#include <set>
#include <string>

#include <lemon/lgf_reader.h>

#include <deregnet/usinglemon.h>
#include <deregnet/utils.h>
#include <deregnet/AvgSbgrphFinder.h>
#include <deregnet/AvgSbgrphModel.h>

using namespace std;

namespace deregnet {

void AvgSbgrphFinder::Data::read_graph(string* pathToLgf) {
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
}

void AvgSbgrphFinder::Data::read_score(string* pathToTsv) {
    if (!graph) {
        cerr << "Error: call AvgSbgrphFinder::Data::read_graph() first." << endl;
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
              score_map[key] = value;
        }
        score = new NodeMap<double>(*graph);
        for (NodeIt v(*graph); v != INVALID; ++v)
            (*score)[v] = score_map[(*nodeid)[v]];
        tsv.close();
    }
    catch ( ... ) {
      cerr << "Could not read score file " << pathToTsv << "." << endl;
      exit(INPUT_ERROR);
    }
}

void AvgSbgrphFinder::Data::get_node_set(std::set<Node>* node_set,
                                         std::set<std::string>* node_ids) {
    if (!graph) {
        cerr << "Error: call AvgSbgrphFinder::Data::read_graph() first." << endl;
        exit(INPUT_ERROR);
    }
    if (node_ids) {
        if (!node_set)
            node_set = new set<Node>();
        Node next;
        for (auto id : *node_ids)
            if ( getNodeById(graph, nodeid, id, &next) )
                node_set->insert(next);
    }
}

void AvgSbgrphFinder::Data::get_root(std::string* root_id) {
    if (root_id) {
        root = new Node();
        if (!getNodeById(graph, nodeid, *root_id, root)) {
            cerr << "Root seems not to be in the graph." << endl;
            exit(INPUT_ERROR);
        }
    }
}

AvgSbgrphFinder::AvgSbgrphFinder(AvgSbgrphFinder::Data xdata,
                                 bool xreceptor_as_root)
                                : data { xdata },
                                  receptor_as_root { xreceptor_as_root }
                               { }

vector<Sbgrph> AvgSbgrphFinder::run(grbfrc::Algorithm algorithm, bool start_heuristic, string model_sense) {
    // set up model
    AvgSbgrphModel model(data.graph, data.score, data.nodeid, data.root);
    model.createVariables();
    model.addBaseConstraints(data.min_size, data.max_size);
    if (data.include)
        model.addIncludeConstraints(data.include);
    if (data.exclude)
        model.addExcludeConstraints(data.exclude);
    if (!receptor_as_root) {
        if (data.receptors)
            model.addReceptorConstraints(data.receptors);
        if (data.terminals)
            model.addTerminalConstraints(data.terminals);

    }
    else if (receptor_as_root) {
        model.reverseGraph();
        if (data.terminals && !data.root)
            model.addReceptorConstraints(data.terminals);
        if (data.receptors)
            model.addTerminalConstraints(data.receptors);
    }

    vector<Sbgrph> subgraphs;
    // find optimal subgraph (modulo time_limit and/or gap_cut)

    if ( model.solve(algorithm,
                     start_heuristic,
                     data.time_limit,
                     data.gap_cut,
                     model_sense) ) {
        subgraphs.push_back ( toSbgrph( model.getCurrentSolution(), "optimal" ) );
    }
    else {
        cerr << "No optimal subgraph could be found." << endl;
        exit(NO_OPTIMAL_SUBGRAPH_FOUND);
    }
    // find suboptimal subgraphs
    set<string> nodes_so_far { subgraphs.back().nodes }; // rethink this policy
    for (int i = 0; i < data.num_subopt_iter; i++) {
        model.addSuboptimalityConstraint(nodes_so_far, data.max_overlap);
        if ( model.solve(algorithm,
                         start_heuristic,
                         data.time_limit,
                         data.gap_cut,
                         model_sense) ) {
            AvgSbgrphModel::Solution current { model.getCurrentSolution() };
            subgraphs.push_back(toSbgrph( current, "suboptimal_" + to_string(i) ));
            for (auto node: current.nodes)
                nodes_so_far.insert(node);
        }
        else {
            cout << "With these parameters, no more than " << i
                 << " suboptimal subgraphs could be found." << endl;
            break;
        }
    }
    return subgraphs;
}

Sbgrph AvgSbgrphFinder::toSbgrph( AvgSbgrphModel::Solution solution, string signature ) {
    Sbgrph subgraph;
    subgraph.signature = signature;
    subgraph.rootid = solution.rootid;
    subgraph.nodes = solution.nodes;
    subgraph.total_score = solution.total_score;
    subgraph.avg_score = solution.avg_score;
    if (data.receptors)
        for (auto nodeid : subgraph.nodes) {
            Node node;
            if (getNodeById(data.graph, data.nodeid, nodeid, &node))
                if (data.receptors->find(node) != data.receptors->end())
                    subgraph.receptors.insert(nodeid);
        }
    if (data.terminals)
        for (auto nodeid : subgraph.nodes) {
            Node node;
            if (getNodeById(data.graph, data.nodeid, nodeid, &node))
                if (data.terminals->find(node) != data.terminals->end())
                    subgraph.receptors.insert(nodeid);
        }
    for (auto edge : solution.edgelist) {
        if (data.receptor_as_root)
            subgraph.edgelist.insert(std::make_pair(edge.second, edge.first));
        else
            subgraph.edgelist = solution.edgelist;
    }
    return subgraph;
}

}   //  namespace deregnet
