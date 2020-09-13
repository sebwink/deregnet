# DeRegNet - Find deregulated subnetworks

![Docker](https://github.com/sebwink/deregnet/workflows/Docker/badge.svg?branch=master)

by [Sebastian Winkler](https://github.com/sebwink) and [Applied Bioinformatics Group](https://kohlbacherlab.org/), University of Tuebingen

## Introduction

One of the main challenges of high-throuput omics technologies 
(genomics, transcriptomics, proteomics, metabolomics, etc.) is the 
interpretation and analysis of the resulting datasets in terms
of known or previously unknown biologial processes. Biological networks
(transcriptional regulatory networks, signaling networks, metabolic
network, etc.) provide promising scaffolds for approaching
multi-omics datasets. Existing resources, constructed for
example from pathway databases like KEGG, Reactome, etc., provide
extensive interconnected networks linking genes, proteins and other
biological agents by various kinds of interactions like generic
activation or inhibition, transcriptional suppression or postranscriptional
modifications like posphorylation. DeRegNet allows the extraction and
prioritisation of subnetworks of larger biomolecular networks based on
suitable omics data like for example gene expression.

## Run via Docker 

Using deregnet via Docker is the only officially supported and documented way of running deregnet. See [examples here](https://github.com/sebwink/deregnet/tree/master/examples).

In any case, you need [Docker](https://www.docker.com/) installed. The next thing you need is a [Gurobi](https://www.gurobi.com/de/) [license](https://www.gurobi.com/de/academia/academic-program-and-licenses/). You can run deregnet with either
a [token server/floating](https://www.gurobi.com/documentation/9.0/quickstart_linux/creating_a_token_server_cl.html#subsection:clientlicensetoken) or [named user](https://www.gurobi.com/documentation/9.0/quickstart_linux/creating_a_new_academic_li.html#subsection:createacademiclicense) license.

### Gurobi floating license

In case of a Gurobi [token server / floating license](https://www.gurobi.com/documentation/9.0/quickstart_linux/creating_a_token_server_cl.html#subsection:clientlicensetoken) you need to make your license file known to the deregnet 
Docker container. Do this by

```sh
export GUROBI_LICENSE=<path to your license file>
```
before running deregnet.

By default the license file will be expected in *~/.licenses/gurobi*.
You need to make sure that the license server configured in your license is reachable from Docker containers
running on your host.

Once the license is configured, the best way to run deregnet is via [*docker/token-server/run*](https://github.com/sebwink/deregnet/tree/master/docker/token-server) script:

```sh
docker/token-server/run <DEREGNET_IMAGE> <CMD>
```
See below for further information about the \<DEREGNET_IMAGE\> and \<CMD\> placeholders.

### Setup for use with a Gurobi named-user license

In case of a Gurobi [named user academic license](https://www.gurobi.com/documentation/9.0/quickstart_linux/creating_a_new_academic_li.html#subsection:createacademiclicense) you also need to make your license file known to the deregnet 
Docker container. Do this by

```sh
export GUROBI_LICENSE=<path to your license file>
```
before running deregnet.

By default the license file will be expected in *~/.licenses/gurobi*.

In order to make a named user license work for deregnet, one additional step is to find the MAC address with respect to which your license
is registered. Do the following before running deregnet:

```
export MAC_ADDRESS_FOR_GUROBI_DOCKER=<YOUR-MAC-ADDRESS>
```

Finding your right \<YOUR-MAC-ADDRESS\> is system-specific, in case of doubt, try all MAC addresses listed by *ifconfig -a* and proceed by trial and error until
your license is accepted while running deregnet (see below).

Once the license is configured, the best way to run deregnet is via [*docker/named-user/run*](https://github.com/sebwink/deregnet/tree/master/docker/named-user) script:

```sh
docker/named-user/run <DEREGNET_IMAGE> <CMD>
```
See below for further information about the \<DEREGNET_IMAGE\> and \<CMD\> placeholders.

### \<DEREGNET_IMAGE\>

Deregnet Docker images are available from [Docker Hub](https://hub.docker.com/repository/docker/sebwink/deregnet) and [GitHub Packages](https://github.com/sebwink/deregnet/packages). Usually, you should be able to just run:

```sh
docker/token-server/run sebwink/deregnet:latest <CMD>
```

To run a specific release of deregnet run for example:


```sh
docker/token-server/run sebwink/deregnet:0.99.999 <CMD>
```

The best way to run with a specific [supported](https://github.com/sebwink/gurobi-docker/blob/master/GUROBI_VERSIONS) Gurobi version is for example like so:


```sh
docker/token-server/run sebwink/deregnet-grb9.0.2:0.99.999 <CMD>
```

```sh
docker/token-server/run sebwink/deregnet-grb8.1.1:latest <CMD>
```

### \<CMD\>
  
deregnet Docker images support multiple commands. The most straight-forward one is to [use](https://github.com/sebwink/deregnet/tree/master/examples/command-line) the main script for deregnet:

```sh
docker/named-user/run sebwink/deregnet:latest avgdrgnt.py --help
```

```sh
usage: avgdrgnt.py [-h] [--include-file INCLUDE_FILE]
                   [--include-genesets INCLUDE_GENESETS] [--include INCLUDE]
                   [--include-id-type INCLUDE_ID_TYPE]
                   [--exclude-file EXCLUDE_FILE]
                   [--exclude-genesets EXCLUDE_GENESETS] [--exclude EXCLUDE]
                   [--exclude-id-type EXCLUDE_ID_TYPE] [--debug]
                   [--absolute-values] --graph GRAPH --scores SCORE_FILE
                   [--default-score DEFAULT_SCORE] [--score-column SCORE_COL]
                   [--score-file-without-header] [--id-column ID_COL]
                   [--sep SEP] [--biomap-mapper ID_MAPPER]
                   [--score-id-type SCORE_ID_TYPE]
                   [--graph-id-type GRAPH_ID_TYPE]
                   [--graph-id-attr GRAPH_ID_ATTR] [--suboptimal SUBOPTIMAL]
                   [--max-overlap-percentage MAX_OVERLAP] [--gap-cut GAP_CUT]
                   [--time-limit TIME_LIMIT] [--model_sense {min,max}]
                   [--output-path OUTPUT] [--flip-orientation]
                   [--min-size MIN_SIZE] [--max-size MAX_SIZE]
                   [--min-num-terminals MIN_NUM_TERMINALS]
                   [--algorithm {GeneralizedCharnesCooper,Dinkelbach,ObjectiveVariableTransform}]
                   [--receptor-file RECEPTOR_FILE]
                   [--receptor-genesets RECEPTOR_GENESETS]
                   [--receptor RECEPTOR] [--receptor-id-type RECEPTOR_ID_TYPE]
                   [--terminal-file TERMINAL_FILE]
                   [--terminal-genesets TERMINAL_GENESETS]
                   [--terminal TERMINAL] [--terminal-id-type TERMINAL_ID_TYPE]

optional arguments:
  -h, --help            show this help message and exit
  --include-file INCLUDE_FILE
                        Path to GMT or GRP file containing genes defining the
                        include layer.
  --include-genesets INCLUDE_GENESETS
                        Comma seperated list of geneset names for include
                        layer,only applicable if GMT file provided.
  --include INCLUDE     Comma seperated list of IDs defining the include
                        layer.
  --include-id-type INCLUDE_ID_TYPE
                        Id-type for include layer genesets. Options: all
                        supported by chosen biomap mapper
  --exclude-file EXCLUDE_FILE
                        Path to GMT or GRP file containing genes defining the
                        exclude layer.
  --exclude-genesets EXCLUDE_GENESETS
                        Comma seperated list of geneset names for exclude
                        layer,only applicable if GMT file provided.
  --exclude EXCLUDE     Comma seperated list of IDs defining the exclude
                        layer.
  --exclude-id-type EXCLUDE_ID_TYPE
                        Id-type for exclude layer genesets. Options: all
                        supported by chosen biomap mapper
  --debug               Debug underlying C++ code with gdb.
  --absolute-values     Whether to take absolute values of the scores.
  --graph GRAPH         A graphml file containing the graph you want to run
                        drgnt with.
  --scores SCORE_FILE   A text file containing the scores. See further options
                        below.
  --default-score DEFAULT_SCORE
                        The score of nodes in the graph which are not scored
                        in your score file. Default: 0.0
  --score-column SCORE_COL
                        Column name of (gene) id in your score file. Default:
                        score
  --score-file-without-header
                        Flag to indicate whether the score file has a header
                        or not.
  --id-column ID_COL    Column name of (gene) id in your score file. Default:
                        id
  --sep SEP             The column seperator in your score file.Options:
                        comma, tab. Default: \t
  --biomap-mapper ID_MAPPER
                        biomap mapper you want to use for id mapping. Default:
                        hgnc
  --score-id-type SCORE_ID_TYPE
                        Which id type do you have in your score file? Options:
                        all thosesupported by the biomap mapper you chose or
                        unspecified. Default: same as graph id type
  --graph-id-type GRAPH_ID_TYPE
                        Which id type does the graph have? Options: all those
                        supportedby the biomap mapper you chose or
                        unspecified. Default: unspecifed i.e. None
  --graph-id-attr GRAPH_ID_ATTR
                        Node attribute which contains the relevant id in the
                        graphml. Default: name
  --suboptimal SUBOPTIMAL
                        Number of suboptimal subgraphs you want to find.
                        (Increases runtime)
  --max-overlap-percentage MAX_OVERLAP
                        How much can suboptimal subgraphs overlap with already
                        found subgraphs. Default: 0
  --gap-cut GAP_CUT     Stop optimization prematurely if current solution
                        within GAP of optimal solution. Default: None
  --time-limit TIME_LIMIT
                        Set a time limit in seconds. Default: None
  --model_sense {min,max}
                        Model sense. Default: max
  --output-path OUTPUT  Folder to which output is written. (Does not have to
                        exist.) Default : cwd
  --flip-orientation    Set --flip-orientation when you want to flip the
                        orientation of the underlying graph.
  --min-size MIN_SIZE   Minimal size of the resulting subgraph(s). Default :
                        15
  --max-size MAX_SIZE   Maximal size of the resulting subgraph(s). Default :
                        15
  --min-num-terminals MIN_NUM_TERMINALS
                        Minimum number of terminals in the resulting
                        subgraph(s). Default : 0
  --algorithm {GeneralizedCharnesCooper,Dinkelbach,ObjectiveVariableTransform}
                        Algorithm to use to solve the fractional integer
                        programming problem.Default: GeneralizedCharnesCooper.
  --receptor-file RECEPTOR_FILE
                        Path to GMT or GRP file containing genes defining the
                        receptor layer.
  --receptor-genesets RECEPTOR_GENESETS
                        Comma seperated list of geneset names for receptor
                        layer,only applicable if GMT file provided.
  --receptor RECEPTOR   Comma seperated list of IDs defining the receptor
                        layer.
  --receptor-id-type RECEPTOR_ID_TYPE
                        Id-type for receptor layer genesets. Options: all
                        supported by chosen biomap mapper
  --terminal-file TERMINAL_FILE
                        Path to GMT or GRP file containing genes defining the
                        terminal layer.
  --terminal-genesets TERMINAL_GENESETS
                        Comma seperated list of geneset names for terminal
                        layer,only applicable if GMT file provided.
  --terminal TERMINAL   Comma seperated list of IDs defining the terminal
                        layer.
  --terminal-id-type TERMINAL_ID_TYPE
                        Id-type for terminal layer genesets. Options: all
                        supported by chosen biomap mapper
```

For example:

```sh 
docker/name-user/run sebwink/deregnet:latest avgdrgnt.py \
  --graph test/kegg_hsa.graphml \
  --scores test/data/score.csv \
  --sep , \
  --graph-id-attr ensembl
``` 

Other commands include *drgnt.py* (Optimization for absolute, not average, best subgraphs).

The most frequent other use cases are to run [Jupyter Lab](https://github.com/sebwink/deregnet/tree/master/examples/jupyterlab) or custom [Python scripts](https://github.com/sebwink/deregnet/tree/master/examples/custom-python-script).

### [Examples](https://github.com/sebwink/deregnet/tree/master/examples)

### Input-output via Docker volumes

Generally your current working directory will be mounted in the running Docker containers. Also some Docker-necessitated access right and owner
sanitations will be carried out, [see for example](https://github.com/sebwink/deregnet/blob/master/docker/named-user/run).

### Security considerations

Only run the deregnet images in trusted environments.

## General remarks

Feedback and problems can be reported via GitHub [Issues](https://github.com/sebwink/deregnet/issues).

[ProjectPage](https://sebwink.github.io/deregnet/)
