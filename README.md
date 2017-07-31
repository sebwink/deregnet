# Deregnet - Find deregulated subnetworks

## Introduction

blup

## Installation

### Prerequisites

This software work is developed and tested on Ubuntu and Scientific Linux, but should
work on other Linux installations as well.

Nothing is offcially supported or guaranteed. Use at your own peril.
Inofficially you can always write a mail to: winkler@informatik.uni-tuebingen.de

You will need two prerequistes installed on your system:

* [Lemon graph library](http://lemon.cs.elte.hu/trac/lemon)
* [Gurobi](http://www.gurobi.com)

### Install

In common.mak set LEMON\_HOME according to your environment.

In gurobi\_version.mak specify your Gurobi version and set GUROBI\_HOME according to your environment.

Then at your command prompt, do the following: 

'''
>> make
'''

## Documentation

[ProjectPage](https://sebwink.github.io/deregnet/)
