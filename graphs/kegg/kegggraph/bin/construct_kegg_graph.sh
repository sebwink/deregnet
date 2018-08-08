#!/bin/bash

PYTHON=python3
RSCRIPT=Rscript
BINDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
SPECIES=$1
KGMLDIR=$BINDIR/../../kgml/$SPECIES
DESTINATION=$BINDIR/../$SPECIES
mkdir -p $DESTINATION
SIFFILE=$DESTINATION/kegg_$SPECIES.sif
if [ ! -d "$KGMLDIR" ]; then
    mkdir $KGMLDIR
    $PYTHON $BINDIR/../../bin/download_kgml.py $SPECIES $KGMLDIR
fi
$RSCRIPT --vanilla $BINDIR/parse_to_graph.R $KGMLDIR $SIFFILE
$PYTHON $BINDIR/to_graphml.py $SIFFILE
