#!/bin/bash

PYTHON=python3
BINDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
SPECIES=$1
KGMLDIR=$BINDIR/../../kgml/$SPECIES
DESTINATION=$BINDIR/../$SPECIES
mkdir -p $DESTINATION
if [ ! -d "$KGMLDIR" ]; then
    mkdir -p $KGMLDIR
    $PYTHON $BINDIR/../../bin/download_kgml.py $SPECIES $KGMLDIR
fi

# GraphML does not work for KEGGtranslators batch mode ...
#$BINDIR/keggtranslator --input $KGMLDIR \
#	               --output $DESTINATION \
#		       --format GraphML \
#		       --remove-orphans

#for kgml in $KGMLDIR/*; do
#    echo $kgml
#    filename=$(basename -- "$kgml")
#    filename="${filename%.*}"
#    $BINDIR/keggtranslator --input $kgml \
#	                   --output $DESTINATION/$filename.graphml \
#			   --format GraphML \
#			   --remove-orphans
#done

for graphml in $DESTINATION/*; do
    echo $graphml
    $PYTHON $BINDIR/make_graphml_igraph_readable.py $graphml
done
#$PYTHON $BINDIR/merge_graphmls.py $DESTINATION
