# Overview



This page is devoted to document and make easier the use of RAFF-
**R**obust **A**lgebraic **F**itting **F**unction. Our intent is to
provide a package to determine fitting functions for a dataset with
ability to detect possible outliers of the dataset. All the code was
made in [Julia language](https://julialang.org), version 1.0.

This package **is not** an implentation of classical least square
solvers. It is an optimization-based package, based on algorithms for Lower Order-Value Optimization (LOVO) which were introduced in [1] and revisited in [2] to fit the user-provided models to experimental data. Recently, a good review can be found in [3]. To find possible outliers, LOVO methods depends on the number of outliers as input information. `RAFF` differs in this point and has no dependence on this number of outliers to perform the fitting process. In order to find a robust adjustment, a voting system is used, which is also responsible for the detection of possible outliers.

## Current Status

The current status of this project is beta quality, don't use for
anything important.  We provide support to serial and parallel
running.

## Developed by

This project was developed by the optimization group at Department of
Mathematics, State University of Maringá, Brazil

* Francisco Sobral (Leader)
* Emerson Vitor Castelani
* Ronaldo Lopes
* Wesley Shirabayashi

## References

[1]Andreani, R., Dunder, C. & Martínez, J.M. Math Meth Oper Res (2005) 61: 365. https://doi.org/10.1007/s001860400410

[2] Andreani, R., Martínez, J. M., Martínez, L., & Yano, F. S. (2009). Low order-value  optimization and applications. *Journal of Global Optimization*, 43(1), 1-22.

[3] Martínez, J.M. TOP (2012) 20:75.https://doi.org/10.1007/s11750-010-0169-1


