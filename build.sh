#!/bin/bash

mkdir -p build
#grbfrc
cd grbfrc && make && cd ..
# drgnt
make -f drgnt.mak
# avgdrgnt
make -f avgdrgnt.mak
# libderegnet
#make -f libderegnet.mak
