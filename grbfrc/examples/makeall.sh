#!/bin/bash

exnum=1

while [ $exnum -le 10 ]; do
  make EX=example${exnum}
  exnum=$((exnum + 1))
done

make clean
