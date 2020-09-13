// --------------------------------------------------------------------------
//                   deregnet -- find deregulated pathways
// --------------------------------------------------------------------------
// Copyright Sebastian Winkler --- Eberhard Karls University Tuebingen, 2020
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
// src/bin/avgdrgnt.cpp
//
// libc #####################################################################
#include <unistd.h>
#include <getopt.h>
#include <sys/types.h>
#include <sys/stat.h>
// libstdc++ libc ###########################################################
#include <cstdlib>
#include <cstring>
// libstdc++ ################################################################
#include <iostream>
#include <string>
#include <set>
// grbfrc ###################################################################
#include "grbfrc/FMILP.hpp"
// deregnet #################################################################
#include "version.hpp"
#include "utils.hpp"
#include "AvgdrgntData.hpp"
#include "DeregnetFinder.hpp"

#define OPTION_ERROR 1

// ##########################################################################

using namespace std;
using namespace grbfrc;
using namespace deregnet;

// Options ##################################################################

struct Options {
	string* lgf_file { nullptr };                      /*< Path to lgf file containing graph to be used */            
    string* score_tsv { nullptr };                     /*< Path to tsv with scores */
    string* root_id { nullptr };                       /*< root node id */
    set<string>* terminals { nullptr };                /*< terminal node ids */
	set<string>* receptors { nullptr };                /*< receptor node ids */
	set<string>* include { nullptr };                  /*< node ids of nodes to be included in subgraph */
	set<string>* exclude { nullptr };                  /*< node ids of nodes to be excluded in subgraph */
	bool no_max_size { false };                        /*< whether to restrict to subgraphs with a maximal size */
    string* outdir { nullptr };                        /*< directory to which to write output and logs */
    bool absolute_values { false };
};

// Data #####################################################################

using Data = AvgdrgntData;

// parse_options ############################################################

void parse_options(int argc, char* argv[], Options& options, Data& data);

// finalize_data ############################################################

void finalize_data(Options& options, Data& data);

// writeSubgraphs ###########################################################

void writeSubgraphs(vector<Subgraph>& subgraphs, string* outdir);

// print_version ############################################################

void print_version();

// print_help ###############################################################

void print_help();

// main #####################################################################

int main(int argc, char* argv[]) {

    Options options;
    Data data;
    parse_options(argc, argv, options, data);
    finalize_data(options, data);
    DeregnetFinder<FMILP, Data> subgraphFinder(&data);
    vector<Subgraph> subgraphs { subgraphFinder.run() };
    std::cout << subgraphs.size() << std::endl;
    writeSubgraphs(subgraphs, options.outdir);
    return 0;
}

// parse_options ############################################################

void parse_options(int argc, char* argv[], Options& options, Data& data) {

	static struct option long_options[] {
		{"graph", 1, 0, 'g'},
		{"terminals", 1, 0, 'T'},
		{"terminals-file", 1, 0, 0},
        {"min-num-terminals", 1, 0, 0},
		{"receptors", 1, 0, 'R'},
		{"receptors-file", 1, 0, 0},
		{"score", 1, 0, 's'},
		{"algorithm", 1, 0, 'a'},
		{"root", 1, 0, 'r'},
		{"flip-orientation", 0, 0, 'f'},
		{"output-dir", 1, 0, 'o'},
		{"min-size", 1, 0, 'l'},
		{"max-size", 1, 0, 'u'},
        {"no-max-size", 0, 0, 0},
		{"time-limit", 1, 0, 't'},
		{"gap-cut", 1, 0, 'c'},
		{"exclude", 1, 0, 'e'},
		{"include", 1, 0, 'i'},
		{"exclude-file", 1, 0, 0},
		{"include-file", 1, 0, 0},
		{"suboptimal", 1, 0, 'b'},
		{"max-overlap-percentage", 1, 0, 'p'},
		{"no-start-heuristic", 0, 0, 0},
        {"help", 0, 0, 'h'},
		{"version", 0, 0, 'v'},
        {"model-sense", 1, 0, 0},
        {"absolute-values", 0, 0, 0},
		{NULL, 0, NULL, 0}
	};

	int opt;
	int option_index = 0;
    string option_signature = "hvg:s:a:l:u:r:R:T:i:e:b:p:t:c:o:";
	while ( (opt = getopt_long(argc, argv, option_signature.c_str(),
	                           long_options, &option_index)) != -1) {
        switch (opt) {
            case 0: {
                string no_start_heuristic { "no-start-heuristic" };
                string no_max_size { "no-max-size" };
                string terminals_file { "terminals-file" };
                string receptors_file { "receptors-file" };
                string include_file { "include-file" };
                string exclude_file { "exclude-file" };
                string model_sense { "model-sense" };
                string absolute_values { "absolute-values" };
                string min_num_terminals { "min-num-terminals" };
                // --no-start-heuristic
                if ( strcmp(long_options[option_index].name, no_start_heuristic.c_str()) == 0 )
                    data.start_heuristic = false;
                // --no-max-size
                else if ( strcmp(long_options[option_index].name, no_max_size.c_str()) == 0)
                    options.no_max_size = true;
                // --terminals-file
                else if ( strcmp(long_options[option_index].name, terminals_file.c_str()) == 0)
                    register_node_set_file(&options.terminals, optarg);
                // --absolute-values
                else if ( strcmp(long_options[option_index].name, absolute_values.c_str()) == 0 )
                    options.absolute_values = true;
                // --receptors-file
                else if ( strcmp(long_options[option_index].name, receptors_file.c_str()) == 0)
                    register_node_set_file(&options.receptors, optarg);
                // --include-file
                else if ( strcmp(long_options[option_index].name, include_file.c_str()) == 0)
                    register_node_set_file(&options.include, optarg);
                // --exclude-file
				else if ( strcmp(long_options[option_index].name, exclude_file.c_str()) == 0)
					register_node_set_file(&options.exclude, optarg);
                else if ( strcmp(long_options[option_index].name, model_sense.c_str()) == 0) 
				    data.model_sense = optarg;
                else if ( strcmp(long_options[option_index].name, min_num_terminals.c_str()) == 0)
                    data.min_num_terminals = new int( atoi(optarg) );
                break;
			}
			case 'g':
				// -g,--graph
				options.lgf_file = new string(optarg); 
				break;
			case 'T':
				// -T,--terminals
				register_node_set(&options.terminals, optarg);
				break;
			case 'R':
				// -R,--receptors
				register_node_set(&options.receptors, optarg);
				break;
			case 's':
				// -s,--score
				options.score_tsv = new string(optarg); 
				break;
			case 'a': {
                // -a,--algorithm
                string gcc { "gcc" };
				string dta { "dta" };
				string ovt { "ovt" };
				if ( strcmp(optarg, gcc.c_str()) == 0 ) 
                    data.algorithm = Algorithm::GCC;
				else if ( strcmp(optarg, dta.c_str()) == 0 ) 
                    data.algorithm = Algorithm::DTA;
				else if ( strcmp(optarg, ovt.c_str()) == 0 ) 
                    data.algorithm = Algorithm::OVT;
				else {
				    cout << "\nUnknown algorithm option provided:" << endl;
					cout << "-a,--algorithm  < gcc | dta | ovt >" << endl;
					exit(OPTION_ERROR);
				}
				break;
			}
			case 'r':
				// -r,--root
				options.root_id = new string(optarg);   
				break;
			case 'f':
				// -f,--flip-orientation
                data.receptor_as_root = false;
				break;
			case 'o':
				// -o,--output-dir
				options.outdir = new string(optarg);
				break;
			case 'l':
				// -l,--min-size
                data.min_size = atoi( optarg );
				break;
			case 'u':
				// -u,--max-size
                data.max_size = atoi( optarg );
				break;
			case 't':
				// -t,--time-limit
				data.time_limit = new double( atof( optarg ) );
				break;
			case 'c':
				// -c,--gap-cut
				data.gap_cut = new double( atof( optarg ) );
				break;
			case 'e':
				// -e,--exclude
				register_node_set(&options.exclude, optarg);
				break;
			case 'i':
				// -i,--include
				register_node_set(&options.include, optarg);
				break;
			case 'b':
				// -b,--suboptimal
				data.num_subopt_iter = atoi( optarg );
				break;
			case 'p':
				// -p,--max-overlap-percentage
				data.max_overlap = atof( optarg );
				break;
			case 'h': 
				// -h,--help
				print_help();
                exit(EXIT_SUCCESS);
			case 'v':
				// -v,--version
                print_version();
                exit(EXIT_SUCCESS);
			case '?':
				break;
        }
	}
}

// finalize_data ############################################################

void finalize_data(Options& options, Data& data) {
    data.read_graph(options.lgf_file);
    data.read_score(options.score_tsv, options.absolute_values);
    data.get_node_set(&data.terminals, options.terminals);
    data.get_node_set(&data.receptors, options.receptors);
    data.get_node_set(&data.include, options.include);
    data.get_node_set(&data.exclude, options.exclude);
    data.get_root(options.root_id);
    if (options.no_max_size)
        data.max_size = countNodes(*data.graph);
    if (data.algorithm == Algorithm::DTA)
        data.gap_cut = nullptr;
}

// writeSubgraphs ###########################################################

void writeSubgraphs(vector<Subgraph>& subgraphs, string* outdirp) {
    // hard to do platform-independent without Boost or C++17
    // on the other hand, it is 2017 now ...
    std::cout << (*outdirp) << std::endl;
    char cwd[2048];
    string outdir;
    if (getcwd(cwd, sizeof(cwd)) != NULL)
        outdir = string(cwd);
    if (outdirp)
        outdir = *outdirp;
    struct stat st;
    if (stat(outdir.c_str(), &st) == -1)
        mkdir(outdir.c_str(), 0700);
    string outdir_plain { outdir + "/plain" };
    if (stat(outdir_plain.c_str(), &st) == -1)
        mkdir(outdir_plain.c_str(), 0700);
    for (auto subgraph : subgraphs) {
        subgraph.writeToFile(outdir);
        subgraph.writeToFile(outdir_plain, true);
    }
}

/*
void write_subgraphs(vector<Subgraph>& subgraphs, string* outdirp) {
    //fs::path outdir(*outdirp);
}
*/

// print_version ############################################################

void print_version() {
    cout << "____________________________________________________\n"
         << "\n|||||||||  avgdrgnt, deregnet version "
         << DEREGNET_VERSION_MAJOR << "."
         << DEREGNET_VERSION_MINOR
         << "  |||||||||\n"
         << "____________________________________________________\n"
         << endl;
}

// print help helpers #######################################################

void print_help_a() {
    cout << " -a,--algorithm:\n\n"
         << " Specify the algorithm with which to solve\n"
         << " the Fractional Integer Linear Program (FILP).\n\n"
         << " Default: gcc.\n\n"
         << " Usage: -a,--algorithm  gcc | dta | ovt \n"
         << " -----------------------------------------\n"
         << " gcc: Generalized Charnes-Cooper transform\n"
         << " dta: Dinkelbach-type algorithm\n"
         << " ovt: Objective-variable transform\n"
         << " ------------------------------------------\n"
         << "____________________________________________________"
         << endl;
}

void print_help_b() {
    cout << "\n -b,--suboptimal:\n\n"
         << " Specify the number of suboptimal subgraphs\n"
         << " to be searched for. \n\n"
         << " Default: 0.\n\n"
         << " Usage: -b,--suboptimal N \n"
         << "         N: integer in [0,100].\n\n"
         << " See also: -p,--max-overlap-percentage\n"
         << "____________________________________________________"
         << endl;
}

void print_help_c() {
    cout << "\n -c,--gap-cut: \n\n"
         << " Specify a gap cut to stop optimization prematurely.\n"
         << " This option is not available with algorithm dta.\n\n"
         << " Default: 0.0 (i.e. solve to optimality).\n\n"
         << " Usage: -c,--gap-cut p\n"
         << "         p: double in [0.0, 100.0].\n\n"
         << " See also: -t,--time-limit\n"
         << "____________________________________________________"
         << endl;
}

void print_help_e() {
    cout << "\n -e,--exclude: \n\n"
         << " Specify a set of nodes which you want to exclude\n"
         << " from any subgraph. Nodes in your set which are not\n"
         << " in the graph are ignored.\n\n"
         << " Default: \"\" (i.e. exclude none of the nodes).\n\n"
         << " Usage: -e,--exclude nodes \n"
         << "         nodes: comma-seperated list of node ids\n"
         << "                (no spaces!).\n\n"
         << " Usage example: -e TP53,DIRAS3,BRCA1 \n\n"
         << " See also: --exclude-file\n"
         << "           -i,--include\n"
         << "           --include-file\n"
         << "____________________________________________________"
         << endl;
}

void print_help_f() {
    cout << "\n -f,--flip-orientation: \n\n"
         << " Apply any algorithm to a reversed version of the\n"
         << " graph. This will result of subgraphs converging\n"
         << " into a terminal node (or several terminal nodes)\n\n"
         << " Usage: -f,--flip-orientation\n"
         << "____________________________________________________"
         << endl;
}

void print_help_g() {
    cout << "\n -g,--graph: \n\n"
         << " Specify the graph with respect to which you want\n"
         << " to look for subgraphs.\n\n"
         << " Default: some version of KEGG\n\n"
         << " Usage: -g,--graph graph.lgf \n"
         << "        graph.lgf: lemon graph format (lgf) file\n"
         << "                   containing your graph. See User's\n"
         << "                   Manual for further specification.\n"
         << "____________________________________________________"
         << endl;
}

void print_help_h() {
    cout << "\n -h,--help: \n\n"
         << " Print this help screen.\n"
         << "____________________________________________________"
         << endl;
}

void print_help_i() {
    cout << "\n -i,--include: \n\n"
         << " Specify a set of nodes which you want to include\n"
         << " in any subgraph. Nodes in your set which are not\n"
         << " in the graph are ignored.\n\n"
         << " Default: \"\" (i.e. no node is included for sure).\n\n"
         << " Usage: -i,--include nodes\n"
         << "         nodes: comma-seperated list of node ids\n"
         << "                (no spaces!).\n\n"
         << " Usage example: -i TP53,DIRAS3,BRCA1 \n\n"
         << " See also: --include-file\n"
         << "           -e,--exclude\n"
         << "           --exclude-file\n"
         << "____________________________________________________"
         << endl;
}

void print_help_l() {
    cout << "\n -l,--min-size: \n\n"
         << " Specify a minimal size of the subgraphs.\n\n"
         << " Default: 20 \n\n"
         << " Usage: -l,--min-size ms \n"
         << "        ms: integer in [1,inf].\n\n"
         << " See also: -u,--max-size\n"
         << "____________________________________________________"
         << endl;
}

void print_help_o() {
    cout << "\n -o,--output-dir: \n\n"
         << " Specify an absolute/relative path of a directory to\n"
         << " which all output will be written.\n\n"
         << " Default: directory from which avgdrgnt is invoked.\n\n"
         << " Usage: -o,--output-dir path/to/directory\n"
         << "____________________________________________________"
         << endl;
}

void print_help_p() {
    cout << "\n -p,--max-overlap-percentage: \n\n"
         << " Specify a percentage which determines to which\n"
         << " extent different suboptimal subgraph can overlap.\n"
         << " TODO: specify exact policy ...\n\n"
         << " Default: 0.0 (i.e. no overlap) \n\n"
         << " Usage: -p,--max-overlap-precentage p \n"
         << "         p: double in [0.0, 100.0].\n\n"
         << " See also: -b,--suboptimal\n"
         << "____________________________________________________"
         << endl;
}

void print_help_r() {
    cout << "\n -r,--root: \n\n"
         << " Fix a particular node to be the root node of any\n"
         << " subgraph. Node identifiers which are not in the \n"
         << " the graph are ignored.\n\n"
         << " Default: None (i.e. the algorithm finds a root)\n\n"
         << " Usage: -r,--root rootid\n"
         << "         rootid: root id\n\n"
         << " Usage example: -r TP53\n\n"
         << " See also: -u,--max-size\n"
         << "____________________________________________________"
         << endl;
}

void print_help_R() {
    cout << "\n -R,--receptors: \n\n"
         << " Specify a set of nodes which you want to be\n"
         << " considered \"receptors\". Nodes in your set\n"
         << " which are not in the graph are ignored.\n\n"
         << " Default: \"\" (i.e. no receptors).\n\n"
         << " Usage: -R,--receptors nodes \n"
         << "         nodes: comma-seperated list of node ids\n"
         << "                (no spaces!).\n\n"
         << " Usage example: -R TP53,DIRAS3,BRCA1 \n\n"
         << " See also: --receptors-file\n"
         << "           -T,--terminals\n"
         << "           --terminals-file\n"
         << "____________________________________________________"
         << endl;
}

void print_help_s() {
    cout << "\n -s,--score: \n\n"
         << " Specify a relative/absolute path to a score file.\n"
         << " Nodes non-existent in the graph are ignored.\n\n"
         << " Default: None (i.e. you have to specify this!)\n\n"
         << " Usage: -s,--score path/to/score_file\n"
         << "         path/to/score_file: score_file should be\n"
         << "         a tsv file (with no header) containing (no-\n"
         << "         thing but) two columns, where the first \n"
         << "         refers to node ids and the second to their\n"
         << "         scores.\n\n"
         << " Use avgdrgnt.py for more flexibility.\n"
         << "____________________________________________________"
         << endl;
}

void print_help_t() {
    cout << "\n -t,--time-limit: \n\n"
         << " Specify a time limit to stop optimization\n"
         << " prematurely.\n\n"
         << " Default:  inf (i.e. no time limit).\n\n"
         << " Usage: -t,--time-limit limit\n"
         << "         limit: double in [0.0, inf].\n\n"
         << " See also: -c,--gap-cut\n"
         << "____________________________________________________"
         << endl;
}

void print_help_T() {
    cout << "\n -T,--terminals: \n\n"
         << " Specify a set of nodes which you want to be\n"
         << " considered \"terminals\". Nodes in your set\n"
         << " which are not in the graph are ignored.\n\n"
         << " Default: \"\" (i.e. no receptors).\n\n"
         << " Usage: -T,--terminals nodes \n"
         << "         nodes: comma-seperated list of node ids\n"
         << "                (no spaces!).\n\n"
         << " Usage example: -T TP53,DIRAS3,BRCA1 \n\n"
         << " See also: --terminals-file\n"
         << "           -R,--receptors\n"
         << "           --receptors-file\n"
         << "____________________________________________________"
         << endl;
}

void print_help_u() {
    cout << "\n -u,--max-size: \n\n"
         << " Specify a maximal size of the subgraphs.\n\n"
         << " Default: 50 \n\n"
         << " Usage: -u,--max-size ms\n"
         << "         ms: integer in [1,inf].\n\n"
         << " See also: -l,--min-size\n"
         << "____________________________________________________"
         << endl;
}

void print_help_v() {
    cout << "\n -v,--version: \n\n"
         << " Display deregnet version information.\n"
         << "____________________________________________________"
         << endl;
}

void print_help_no_start_heuristic() {
    cout << "\n --no-start-heuristic: \n\n"
         << " Do not use greedy start heuristic.\n"
         << "____________________________________________________"
         << endl;
}

void print_help_no_max_size() {
    cout << "\n --no-max-size: \n\n"
         << " Do not restrict the maximal number of nodes in any\n"
         << " subgraph. This option overwrites anything set by \n"
         << " the -u,--max-size option.\n\n"
         << " See also: -u,--max-size\n"
         << "____________________________________________________"
         << endl;
}

void print_help_terminals_file() {
    cout << "\n --terminals-file: \n\n"
         << " Specify a relative/absolute path to a file cont-\n"
         << " aining \"terminals\".\n\n"
         << " Usage: --terminals-file path/to/file\n"
         << "         file: text file containing one node id per\n"
         << "               line (and nothing else).\n\n"
         << " See also: -T,--terminals\n"
         << "           -R,--receptors\n"
         << "           --receptors-file\n"
         << "____________________________________________________"
         << endl;
}

void print_help_receptors_file() {
    cout << "\n --receptors-file: \n\n"
         << " Specify a relative/absolute path to a file cont-\n"
         << " aining \"receptors\".\n\n"
         << " Usage: --receptors-file path/to/file\n"
         << "         file: text file containing one node id per\n"
         << "               line (and nothing else).\n\n"
         << " See also: -R,--receptors\n"
         << "           -T,--terminals\n"
         << "           --terminals-file\n"
         << "____________________________________________________"
         << endl;
}

void print_help_include_file() {
    cout << "\n --include-file: \n\n"
         << " Specify a relative/absolute path to a file cont-\n"
         << " aining nodes to be included in any subgraph.\n\n"
         << " Usage: --include-file path/to/file\n"
         << "         file: text file containing one node id per\n"
         << "               line (and nothing else).\n\n"
         << " See also: -i,--include\n"
         << "           -e,--exclude\n"
         << "           --exclude-file\n"
         << "____________________________________________________"
         << endl;
}

void print_help_exclude_file() {
    cout << "\n --exclude-file: \n\n"
         << " Specify a relative/absolute path to a file cont-\n"
         << " aining nodes to be excluded from any subgraph.\n\n"
         << " Usage: --exclude-file path/to/file\n"
         << "         file: text file containing one node id per\n"
         << "               line (and nothing else).\n\n"
         << " See also: -e,--exclude\n"
         << "           -i,--include\n"
         << "           --include-file\n"
         << "____________________________________________________"
         << endl;
}

// print_help ###############################################################

void print_help() {
    print_version();
    print_help_a();
    print_help_b();
    print_help_c();
    print_help_e();
    print_help_f();
    print_help_g();
    print_help_h();
    print_help_i();
    print_help_l();
    print_help_o();
    print_help_p();
    print_help_r();
    print_help_R();
    print_help_s();
    print_help_t();
    print_help_T();
    print_help_u();
    print_help_v();
    print_help_no_start_heuristic();
    print_help_no_max_size();
    print_help_terminals_file();
    print_help_receptors_file();
    print_help_include_file();
    print_help_exclude_file();
    print_version();
    cout << endl;
}

// ##########################################################################
