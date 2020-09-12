# DeRegNet - Find deregulated subnetworks

![Docker](https://github.com/sebwink/deregnet/workflows/Docker/badge.svg?branch=master)

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

Using deregnet is the only officially supported and documented way of running deregnet.


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
