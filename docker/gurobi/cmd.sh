#!/usr/bin/bash

set -e

export GUROBI_HOME=/opt/gurobi
export LD_LIBRARY_PATH=/opt/gurobi/lib:$LD_LIBRARY_PATH
export PATH=/opt/gurobi/bin:$PATH

exec gurobi.sh
