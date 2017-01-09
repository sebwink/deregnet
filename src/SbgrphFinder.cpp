#include <iostream>
#include <fstream>
#include <set>
#include <string>
#include <map>
#include <utility>

#include <lemon/lgf_reader.h>

#include <deregnet/usinglemon.h>
#include <deregnet/utils.h>
#include <deregnet/SbgrphFinder.h>
#include <deregnet/SbgrphModel.h>

using namespace std;

namespace deregnet {

void SbgrphFinder::Data::read_graph(string* pathToLgf) {
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

void SbgrphFinder::Data::read_score(string* pathToTsv) {
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

void SbgrphFinder::Data::get_node_set(std::set<Node>* node_set,
                                      std::set<std::string>* node_ids) {
    if (!graph) {
        cerr << "Error: call SbgrphFinder::Data::read_graph() first." << endl;
        exit(INPUT_ERROR);
    }
    if (node_ids) {
        if (!node_set)
            node_set = new set<Node>();
        Node next;
        for (auto id : *node_ids) {
            if (!receptor_as_root) {
                if ( getNodeById(revgraph, revnodeid, id, &next) )
                    node_set->insert(next);
            }
            else {
                if ( getNodeById(graph, nodeid, id, &next) )
                    node_set->insert(next);
            }
        }
    }
}

void SbgrphFinder::Data::get_root(std::string* root_id) {
    if (root_id) {
        if (!receptor_as_root) {
            revroot = new Node();
            if (!getNodeById(revgraph, revnodeid, *root_id, revroot)) {
                cerr << "Root seems not to be in the graph." << endl;
                exit(INPUT_ERROR);
            }
        }
        else {
            root = new Node();
            if (!getNodeById(graph, nodeid, *root_id, root)) {
                cerr << "Root seems not to be in the graph." << endl;
                exit(INPUT_ERROR);
            }
        }
    }
}

void SbgrphFinder::Data::get_reversed_graph() {
    receptor_as_root = false;
    map<Node,Node> node_correspondence;
    revgraph = new Graph();
    revscore = new NodeMap<double>(*revgraph);
    revnodeid = new NodeMap<string>(*revgraph);
    for (NodeIt v(*graph); v != INVALID; ++v) {
        Node revv = revgraph->addNode();
        (*revscore)[revv] = (*score)[v];
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
}

SbgrphFinder::SbgrphFinder(SbgrphFinder::Data xdata)
                           : data { xdata }
                          { }

vector<Sbgrph> SbgrphFinder::run(bool start_heuristic, std::string model_sense) {
    // set up model
    SbgrphModel* model;
    if (!data.receptor_as_root)
        model = new SbgrphModel(data.revgraph, data.revscore, data.revnodeid, data.revroot);
    else
        model = new SbgrphModel(data.graph, data.score, data.nodeid, data.root);

    model->createVariables();
    model->addBaseConstraints(data.size);
    if (data.include)
        model->addIncludeConstraints(data.include);
    if (data.exclude)
        model->addExcludeConstraints(data.exclude);
    if (!data.receptor_as_root) {
        if (data.receptors)
            model->addReceptorConstraints(data.receptors);
        if (data.terminals)
            model->addTerminalConstraints(data.terminals);

    }
    else if (data.receptor_as_root) {
        if (data.terminals && !data.root)
            model->addReceptorConstraints(data.terminals);
        if (data.receptors)
            model->addTerminalConstraints(data.receptors);
    }

    vector<Sbgrph> subgraphs;
    // find optimal subgraph (modulo time_limit and/or gap_cut)

    if ( model->solve(start_heuristic,
                      data.time_limit,
                      data.gap_cut,
                      model_sense) ) {
        subgraphs.push_back ( toSbgrph( model->getCurrentSolution(), "optimal" ) );
    }
    else {
        cerr << "No optimal subgraph could be found." << endl;
        exit(NO_OPTIMAL_SUBGRAPH_FOUND);
    }

    // find suboptimal subgraphs
    set<string> nodes_so_far { subgraphs.back().nodes };
    for (int i = 0; i < data.num_subopt_iter; i++) {
        model->addSuboptimalityConstraint(nodes_so_far, data.max_overlap, data.size);
        if ( model->solve(start_heuristic,
                          data.time_limit,
                          data.gap_cut,
                          model_sense) ) {
            SbgrphModel::Solution current { model->getCurrentSolution() };
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

Sbgrph SbgrphFinder::toSbgrph( SbgrphModel::Solution solution, string signature ) {
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
        if (!data.receptor_as_root)
            subgraph.edgelist.insert(std::make_pair(edge.second, edge.first));
        else
            subgraph.edgelist = solution.edgelist;
    }
    return subgraph;
}

}   //  namespace deregnet
