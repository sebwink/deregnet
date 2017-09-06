#!/bin/bash

PYTHON="/opt/anaconda/3.5/bin/python"
RSCRIPT="/home/sebastian/tools/R/3.3.2/bin/Rscript"

species=$1
sif_file="kegg_${species}.sif"

mkdir -p kgml
$PYTHON download_kgml.py $species kgml
$RSCRIPT --vanilla parse_to_graph.R kgml $sif_file
$PYTHON to_graphml.py $sif_file
rm -r kgml
mv kegg_${species}.sif ..
mv kegg_${species}.graphml ..
