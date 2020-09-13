#!/usr/bin/env bash

export GUROBI_LICENSE_MODE=named-user
export DEREGNET_HOME=../..

$DEREGNET_HOME/docker/$GUROBI_LICENSE_MODE/run sebwink/deregnet:latest jupyter lab
