---
title: 'RAFF: Robust Algebraic Fitting Function in Julia'
tags:
  - Julia
  - Statistics
  - Lower Order-Value Optimization
  - Outlier detection
  - Nonlinear optimization
authors:
  - name: Emerson V. Castelani
    orcid: 0000-0000-0000-0000
    affiliation: 1
  - name: Ronaldo Lopes
    orcid: 0000-0000-0000-0000
    affiliation: 1
  - name: Wesley V. I. Shirabayashi
    orcid: 0000-0000-0000-0000
    affiliation: 1
  - name: Francisco N. C. Sobral
    orcid: 0000-0003-4963-0946
    affiliation: 1
affiliations:
 - name: Department of Mathematics, State University of Maringá, Paraná, Brazil
   index: 1
date: 11 March 2019
bibliography: paper.bib
---

# Summary

Let $f : \Re^n \to \Re$ be a function whose mathematical description
is not available. This function can be, for example, a black-box, a
proprietary computer program or an experiment. Suppose that a data set
$S = \{(t_1, f(t_1)), \dots, (t_m, f(t_m))\}$ of $m$ observations of
$f$ is available and we want to approximate $f$ by a known model
$\phi$. Model $\phi$ can be defined as $\phi(x;t)$, where $t$ are the
$n$ independent variables of $f$ and $x$ represents some parameters of
$\phi$. ``RAFF`` (Robust Algebraic Fitting Function) is a Julia
package developed to find the optimal parameters $x$ for $\phi$ in
order to adjust it to the observed values $S$ of the unknown function
$f$.

The adjustment of mathematical functions to data is a problem that
appears in many areas of science. When data comes from real
experiments, non-expected errors may cause the appearance of
outliers. Detection of outliers is always regarded as the statistical
part of data adjustment. ``RAFF`` provides automatic detection of
outliers using a voting system. It is an optimization-based package,
based on algorithms for Lower Order-Value Optimization (LOVO)
[@Andreani2009] to fit the user-provided models $\phi$ to experimental
data. In order to find a robust adjustment, a voting system is used,
which is also responsible for the detection of possible outliers.

# Functionality

``RAFF`` main methods expect as input a data set of the observed data
and a model function, whose parameters one intends to adjust. The
model function is a regular Julia function with 2 arguments: $x$
represents the parameters of the model and $t$ represents the
arguments of function $f$. The following function is an example of a
model representing a circle
$$
\phi(x, t) = (t_1 - x_1)^2 + (t_2 - x_2)^2 - x_3^2.
$$

Note that $x_1$ and $x_2$ are parameters related with the center of
the circle, while $x_3$ is related to its radius. Therefore, this
model contains 3 parameters. Since we think that the unknown function
$f$ is a circle, it is obvious that its arguments are a
two-dimensional vector $t$. The observed data can be represented by
the following table:

-------------------------------------
$t_1$       $t_2$     $f(t_1, t_2)$
-------     -------   ---------------
 6.11129    1.0000    0.0
 6.12134    1.0038    0.0
 5.73863    2.1833    0.0
 5.35241    3.5871    0.0
 4.76209    4.2487    0.0
 3.57131    5.3119    0.0
 2.75896    5.7070    0.0
 2.07888    5.7852    0.0
 1.65623    6.0197    0.0
-0.31324    5.7656    0.0
-2.18214    4.9540    0.0
-2.45642    4.4935    0.0
-3.15965    3.8885    0.0
-3.77514    2.5237    0.0
-4.03242    1.4327    0.0
-3.93137   -0.0512    0.0
-2.55007   -2.6777    0.0
 0.16087   -4.0584    0.0
 1.83171   -3.9190    0.0
 4.19091   -2.7741    0.0
 5.67758   -0.6784    0.0
 5.81633    0.2874    0.0
 5.94573    0.3777    0.0
-1.12132    3.1213    0.0
-3.59619   -3.5961    0.0
-------------------------------------

# Additional features

 The user may also provide more information, such as the
expected number of *trusted* observations. Additional methods and
options are also available to more advanced users, such as generation
of random test data and multistart strategies. It uses Julia's
``ForwardDiff``[@ForwardDiff] package, so the user does not need to
compute the derivatives of the model function.

``RAFF`` can be run in serial, parallel and distributed
environments. Parallel and distributed methods use the native
``Distributed``[@Distributed] package. The distributed version is a
centralized implementation that does not use shared arrays, therefore,
can be run both locally or in a cluster of computers.

This package is intended to be used by all experimental researchers
who know a little about mathematical modeling and fitting functions.

# References
