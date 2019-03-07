---
title: 'RAFF: Robust Algebraic Fitting Function in Julia'
tags:
  - Julia
  - Statistics
  - Lower Order-Value Optimization
  - Outlier detection
  - Nonlinear optimization
authors:
  - name: Francisco N. C. Sobral
    orcid: 0000-0003-4963-0946
    affiliation: 1
  - name: Emerson V. Castelani
    orcid: 0000-0000-0000-0000
    affiliation: 1
affiliations:
 - name: Department of Mathematics, State University of Maringá, Paraná, Brazil
   index: 1
date: 26 February 2019
bibliography: paper.bib
---

# Summary

The adjustment of mathematical functions to data is a problem that
appears in many areas of science. When data comes from real
experiments, non-expected error may cause the appearance of
outliers. Detection of outliers is always regarded as the statistical
part of data adjustment. ``RAFF`` (Robust Algebraic Fitting Function)
is a Julia package to the adjustment of mathematical models with
automatic detection of outliers. It is an optimization-based package,
based on algorithms for Lower Order-Value Optimization (LOVO)
[@Andreani2009] to fit the user-provided models. In order to find a
robust adjustment, a voting system is used. The voting system is also
responsible for the detection of possible outliers.

``RAFF`` main methods expect as input a dataset of the observed data
and a model function, whose parameters one intends to adjust. The user
may also provide more information, such as the expected number of
*trusted* observations. Additional methods and options are also
available to more advanced users, such as generation of random test
data and multistart strategies. It uses Julia's
``ForwardDiff``[@ForwardDiff] package, so the user does not need to
compute the derivatives of the model function. ``RAFF`` also has a
parallel and distributed version, which uses the native
``Distributed``[@Distributed] package. The distributed version is a
centralized implementation that does not uses shared arrays,
therefore, can be run both locally or in a cluster of heterogeneous
computers.

``RAFF`` is intended to be used by all experimental researchers who
know a little about mathematical modeling and fitting functions.

# Mathematics

Single dollars ($) are required for inline mathematics e.g. $f(x) = e^{\pi/x}$

Double dollars make self-standing equations:

$$\Theta(x) = \left\{\begin{array}{l}
0\textrm{ if } x < 0\cr
1\textrm{ else}
\end{array}\right.$$


# Citations

Citations to entries in paper.bib should be in
[rMarkdown](http://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html)
format.

# Figures

Figures can be included like this: ![Example figure.](figure.png)

# Acknowledgements

We acknowledge contributions from Brigitta Sipocz, Syrtis Major, and Semyeong
Oh, and support from Kathryn Johnston during the genesis of this project.

# References
