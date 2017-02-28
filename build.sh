#!/bin/bash

mkdir -p build
# drgnt
make -f drgnt.mak
#grbfrc
cd grbfrc && make -f grbfrc.mak && cd ..
# avgdrgnt
make -f avgdrgnt.mak

