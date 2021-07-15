# Overview



This page is devoted to document and make easier the use of RAFF-
**R**obust **A**lgebraic **F**itting **F**unction. Our intent is to
provide a package to determine fitting functions for a dataset with
ability to detect possible outliers of the dataset. All the code was
made in [Julia language](https://julialang.org), version 1.0.

This package **is not** an implentation of classical least squares
solvers. It is an optimization-based package, based on algorithms for Lower Order-Value Optimization (LOVO) which were introduced in [1] and revisited in [2] to fit the user-provided models to experimental data. Recently, a good review can be found in [3]. To find possible outliers, LOVO methods depend on the number of outliers as input information. `RAFF` differs in this point and has no dependence on this number of outliers to perform the fitting process. In order to find a robust adjustment, a voting system is used, which is also responsible for the detection of possible outliers.

## Current Status

The current status of this project is beta quality, don't use for
anything important.  We provide support to serial and parallel
running.

## Developed by

This project was developed by the optimization group at Department of
Mathematics, State University of Maringá, Brazil.

* Francisco Sobral (Leader)
* Emerson Vitor Castelani
* Ronaldo Lopes
* Wesley Shirabayashi

The authors of this package were sponsored by **Fundação Araucária**,
project number 002/17 - 47223.

## References

[1] Andreani, R., Dunder, C. & Martínez, J.M. Math Meth Oper Res (2005) 61: 365. https://doi.org/10.1007/s001860400410

[2] Andreani, R., Martínez, J.M., Martínez, L. et al. J Glob Optim (2009) 43: 1. https://doi.org/10.1007/s10898-008-9280-3

[3] Martínez, J.M. TOP (2012) 20: 75. [https://doi.org/10.1007/s11750-010-0169-1]

## Citing this package

If you would like to cite this package, please use

> Castelani, E. V., Lopes, R., Shirabayashi, W., & Sobral,
> F. N. C. (2019). RAFF.jl: Robust Algebraic Fitting Function in
> Julia. *Journal of Open Source Software*,
> 4(39), 1385. https://doi.org/10.21105/joss.01385

[BibTex](assets/raff.bib)

The following paper describes the theory and several comparison tests

> @article{Castelani2021,
>    author = {Castelani, Emerson V. and Lopes, Ronaldo and Shirabayashi, Wesley V. I. and Sobral, Francisco N. C.},
>    doi = {10.1007/s10898-020-00970-4},
>    journal = {Journal of Global Optimization},
>    title = {{A robust method based on LOVO functions for solving least squares problems}},
>    url = {http://link.springer.com/10.1007/s10898-020-00970-4},
>    year = {2021}
> }
