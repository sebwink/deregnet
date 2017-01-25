/*! \mainpage libgrbfrc documentation

Introduction
============

Libgrbfrc is a C++ library which allows solving fractional mixed-integer linear programs (FMILPs), i.e. mathematical programs with with rational objective where enumerator and denominator are linear, subject to linear constraints where integer and continuous variables are allowed. Libgrbfrc is build on top of the C++ API of the Gurobi optimizer.

A FMILP is given by an optimization problem of the following form:

\f[

\begin{equation}
	\begin{array}{ll}
		\displaystyle \max \:\: \frac{c^Tx + d}{e^Tx + f} \\\\
				      \textrm{s.t.} 
					                & x\in \mathbb{R}^{n_c} \times \mathbb{Z}^{n_d} \\
						            & Ax \leq b
	\end{array}
\end{equation}

\f]

where \f$c,d \in \mathbb{R}^{n}\f$, \f$e,f \in \mathbb{R}\f$, \f$A \in \mathbb{R}^{m \times n}\f$ and \f$b \in \mathbb{R}^{m}\f$ with \f$m \in \mathbb{N}\f$ being the number of (linear) constraints and \f$n = n_c + n_d\f$ being the total number of variables where \f$n_c \in \mathbb{N}\f$ denotes the number of continuous variables and \f$n_d \in \mathbb{N}\f$ denotes the number of integer variables.

Algorithms
==========

Charnes-Cooper transform
------------------------

... a classical approach for solving continuous problems, i.e. fractional linear programs (FLPs).

bla

Generalized Charnes-Cooper transform
------------------------------------

... a reformulation-linearization method to solve general FMILPs.

Objective-Variable transform
----------------------------

... a reformulation-linearization method to solve FMILPs with no continuous variables in the objective enumerator.

Dinkelbach-type algorithm
-------------------------

... a iterative scheme for solving general FMILPs.


Gurobi
======


Installation
============

Find out yourself

Linux
-----

make

OS X
----

Windows
-------

don't


Examples
========



Applications
============
