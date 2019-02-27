---
title: 'RAFF: Robust Algebraic Fitting Function in Julia'
tags:
  - Julia
  - Statistics
  - Lower Order-Value Optimization
  - Outlier detection
  - Nonlinear optimization
authors:
  - name: Francisco N. C Sobral
    orcid: 0000-0003-4963-0946
    affiliation: 1
  - name: Author 2
    orcid: 0000-0000-0000-0000
    affiliation: 2
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
part of data adjustment. ``RAFF`` (Robust Algebraic Fitting Function) is
a Julia package to the adjustment of mathematical models with
automatic detection of outliers.

``RAFF`` is based on algorithms for Lower Order-Value Optimization
(LOVO) [@Andreani2009] to fit the given models. In order to find a
robust adjustment, a voting system is used. The voting system is also
responsible for the detection of possible outliers.

``Gala`` is an Astropy-affiliated Python package for galactic dynamics. Python
enables wrapping low-level languages (e.g., C) for speed without losing
flexibility or ease-of-use in the user-interface. The API for ``Gala`` was
designed to provide a class-based and user-friendly interface to fast (C or
Cython-optimized) implementations of common operations such as gravitational
potential and force evaluation, orbit integration, dynamical transformations,
and chaos indicators for nonlinear dynamics. ``Gala`` also relies heavily on and
interfaces well with the implementations of physical units and astronomical
coordinate systems in the ``Astropy`` package [@astropy] (``astropy.units`` and
``astropy.coordinates``).

``Gala`` was designed to be used by both astronomical researchers and by
students in courses on gravitational dynamics or astronomy. It has already been
used in a number of scientific publications [@Pearson:2017] and has also been
used in graduate courses on Galactic dynamics to, e.g., provide interactive
visualizations of textbook material [@Binney:2008]. The combination of speed,
design, and support for Astropy functionality in ``Gala`` will enable exciting
scientific explorations of forthcoming data releases from the *Gaia* mission
[@gaia] by students and experts alike. The source code for ``Gala`` has been
archived to Zenodo with the linked DOI: [@zenodo]

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
