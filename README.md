# DeRegNet - Find deregulated subnetworks

## Introduction

One of the main challenges of high-throuput omics technologies 
(genomics, transcriptomic, proteomic, metabolomics, etc.) is the 
interpretation and analysis of the resulting datasets in terms
of known or previously unknown biologial processes. Biological networks
(transcriptional regulatory networks, signaling networks, metabolic
network, etc.) provide promising scafolds with which to approach
multi-omics datasets. Existing resources, constructed for
example from pathway databases like KEGG, Reactome, etc., provide
extensive interconnected networks linking genes, proteins and other
biological agents by various kinds of interactions like generic
activation or inhibition, transcriptional suppression or postranscriptional
modifications like posphorylation. DeRegNet allows the extraction and
prioritisation of subnetworks of larger biomolecular networks based on
suitable omics data like for example gene expression.

## Installation

There are, broadly speaking, two modes of installation and usage for the basic
DeRegNet programs: *Native* and *Docker*. It is recommended to use Docker
whenever possible. Both ways assume a Linux operating system.

### Native

#### Prerequisites

DeRegNet relies on the following two software distributions:

* [Lemon graph library](http://lemon.cs.elte.hu/trac/lemon)
* [Gurobi](http://www.gurobi.com) (version >= 7.0.2)

To install DeRegNet locally on your machine install these two prerequsistes first.
DeRegNet requires a Gurobi version >= 7.0.2. The native local installation of 
DeRegNet works with any [supported licensing scheme](http://www.gurobi.com/downloads/licenses/license-center)
for the Gurobi libraries.

#### Install

In *common.mak* set *LEMON\_HOME* according to your environment.

In *gurobi\_version.mak* specify your Gurobi version and set *GUROBI\_HOME* according to your environment.

```sh
make
```

This builds the main DeRegNet executables in the *bin* directory: *drgnt* and *avgdrgnt*.

In order to install the Python interface of DeRegNet navigate to the *python* subdirectory
and run

```sh
python3 setup.py install
```

#### Basic usage

```sh
bin/avgdrgnt.py --help
```

#### Examples

```sh
bin/avgdrgnt.py --graph test/kegg_hsa.graphml \
                --scores test/data/score.csv \
				--sep , \
				--output-path test
```

### Docker

This mode of installation and usage does not require you to download
this repository.

#### Prerequisites

* [Docker](https://www.docker.com/).
* A floating (!) [license for Gurobi](http://www.gurobi.com/downloads/licenses/license-center).
  (currently the corresponding license server is assumed to be called "license")

#### Install

```sh
docker pull sebwink/deregnet
```

#### Basic usage

```sh
docker run sebwink/deregnet --help
```

#### Examples

```sh
docker run -v $(pwd)/test:/io sebwink/deregnet \
                                  --graph kegg_hsa.graphml \
                                  --scores data/score.csv \
				                  --sep , 
```

## General remarks

Feedback and problems can be reported to deregnet@informatik.uni-tuebingen.de.

## Documentation

```sh
bin/avgdrgnt.py --help
```

or

```sh
docker run sebwink/deregnet --help
```

[ProjectPage](https://sebwink.github.io/deregnet/)
