# Deregnet - Find deregulated subnetworks

## Introduction

blup

## Installation

### Prerequisites

This software work is developed and tested on Ubuntu and Scientific Linux, but should
work on other Linux installations as well.

Nothing is offcially supported or guaranteed. Use at your own peril.

Inofficially you can always write a mail to: *winkler@informatik.uni-tuebingen.de*

You will need two prerequistes installed on your system:

* [Lemon graph library](http://lemon.cs.elte.hu/trac/lemon)
* [Gurobi](http://www.gurobi.com)

### Install

In *common.mak* set *LEMON\_HOME* according to your environment.

In *gurobi\_version.mak* specify your Gurobi version and set *GUROBI\_HOME* according to your environment.
If you use Gurobi version < 7.0.2 write me an E-Mail. Better, install the newest version of Gurobi.

Use *bash* otherwise things may break...

Then at your command prompt, type *make*. 

In case of successful compilation you find the two binaries *drgnt* and *avgdrgnt* in the *bin* directory.

## Basic usage

You can either use *bin/drgnt* or *bin/avgdrgnt* directly, or you can use the Python (Python 3!) wrappers
*bin/drgnt.py* and *bin/avgdrgnt.py*. For the Python scripts to be functional make a soft link
to python/deregnet in the relevant "site-packages" directory of the relevant Python installation (regrettably 
there is no canonical Python installation procedure for the "deregnet" package at this point). Also make sure
you point the shebang in *bin/drgnt.py* and *bin/avgdrgnt.py* to the Python installation in question. Also, make 
sure you have the Python packages referenced in *python/deregnet/deregnet.py* installed.

## General remarks

Good luck!

## Documentation

[ProjectPage](https://sebwink.github.io/deregnet/)
