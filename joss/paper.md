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
    orcid: 0000-0001-9718-6486
    affiliation: 1
  - name: Ronaldo Lopes
    affiliation: 1
  - name: Wesley V. I. Shirabayashi
    orcid: 0000-0002-7790-6703	
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

Let $f : \mathbb{R}^n \to \mathbb{R}$ be a function whose mathematical description
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

Linear and nonlinear regression is essentially the adjustment of
mathematical functions to data and is a problem that appears in many
areas of science. When data comes from real experiments, non-expected
errors may cause the appearance of outliers, which might be
responsible for causing the regression calculated by sum of deviations
to result in misleading approximations. Regression is strongly
connected to Statistics but practical methods to detect outliers are
not very common. In [@Motulsky2006a], for example, the authors develop
a method for outlier detection based on the assumption that the error
follows a Lorentzian distribution around the function and use
nonlinear regression based on least squares. ``RAFF.jl`` provides
automatic detection of outliers using a voting system. It is an
optimization-based package, based on algorithms for Lower Order-Value
Optimization (LOVO) which were introduced in [@Andreani2005] and
revisited in [@Andreani2009] to fit the user-provided models $\phi$ to
experimental data. Recently, a good review can be found in
[@Martinez2012]. In order to find a robust adjustment, a voting system
is used, which is also responsible for the detection of possible
outliers.

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

Note that $x_1$ and $x_2$ are parameters related to the center of the
circle, while $x_3$ is related to its radius. Therefore, this model
contains 3 parameters. Since we think that the unknown function $f$ is
a circle, it is obvious that its arguments are given by a
two-dimensional vector $t$. The observed data can be represented by
the following table:

  $t_1$     $t_2$     $f(t_1, t_2)$
-------   -------   ---------------
 6.1112    1.0000            0.0000
 6.1213    1.0038            0.0000
 5.7386    2.1833            0.0000
 5.3524    3.5871            0.0000
 4.7620    4.2487            0.0000
 3.5713    5.3119            0.0000
 2.7589    5.7070            0.0000
 2.0788    5.7852            0.0000
 1.6562    6.0197            0.0000
-0.3132    5.7656            0.0000
-2.1821    4.9540            0.0000
-2.4564    4.4935            0.0000
-3.1596    3.8885            0.0000
-3.7751    2.5237            0.0000
-4.0324    1.4327            0.0000
-3.9313   -0.0512            0.0000
-2.5500   -2.6777            0.0000
 0.1608   -4.0584            0.0000
 1.8317   -3.9190            0.0000
 4.1909   -2.7741            0.0000
 5.6775   -0.6784            0.0000
 5.8163    0.2874            0.0000
 5.9457    0.3777            0.0000
-1.1213    3.1213            0.0000
-3.5961   -3.5961            0.0000

In this example, we are trying to find a circle to fit the observed
data. Therefore, all the values of $f$ should be zero. The observed
data is shown in [Figure 1](#circle). The dots represent the observed
data, the red ones being the outliers. They were generated as a
perturbation of the points lying in the blue dashed circle. Using the
Least Squares technique with the model above, the red circle is
found. When RAFF is applied to the same problem, it correctly
identifies the two outliers. The resulting circle is depicted as the
green circle, very close the true circle.

![Points representing a circle. The red dots are two outliers that
 should be ignored. The blue dashed circle is the true one, while the
 red was obtained by traditional Least Squares techniques and the
 green one was obtained by RAFF.](circle.png){#circle width=50%,
  height=50%}

# Additional features

The user may also provide more information to ``RAFF``, such as the
expected number of *trusted* observations. Additional methods and
options are also available to more advanced users, such as generation
of random test data and multistart strategies. Derivatives of the
model $\phi$ can also be provided, what results in a faster executing
time. When they are not provided by the user, ``RAFF`` uses Julia's
``ForwardDiff`` [@Revels2016] package.

``RAFF`` can be run in serial, parallel and distributed
environments. Parallel and distributed methods use the native
[``Distributed``](https://docs.julialang.org/en/v1.0/stdlib/Distributed/)
package. The distributed version is a centralized implementation that
does not use shared arrays, therefore, can be run both locally or in a
cluster of computers.

This package is intended to be used by all experimental researchers
who know a little about mathematical modeling and fitting functions.

# Instalation and usage

``RAFF`` is an open-source software that can be [downloaded from
Github](https://github.com/fsobral/RAFF.jl). All the description for
first time usage or its API is available at its
[documentation](https://fsobral.github.io/RAFF.jl/stable/).

# References
