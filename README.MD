RAFF - Robust Algebraic Fitting Function
========================================

[![Build Status](https://travis-ci.com/fsobral/RAFF.jl.svg?token=AFD4xpyZzz3My5T6Sx1k&branch=master)](https://travis-ci.com/fsobral/RAFF.jl)
[![codecov](https://codecov.io/gh/fsobral/RAFF.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/fsobral/RAFF.jl)
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://fsobral.github.io/RAFF.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://fsobral.github.io/RAFF.jl/dev)

[![DOI](https://www.zenodo.org/badge/147533027.svg)](https://www.zenodo.org/badge/latestdoi/147533027)
[![DOI](http://joss.theoj.org/papers/10.21105/joss.01385/status.svg)](https://doi.org/10.21105/joss.01385)

This package implements a robust method to fit a given function
(described by its parameters) to input data. The method is based on
the LOVO algorithm [1] and also in a suitable voting strategy in order
automatically eliminate outliers.

Potential use of this package is in statistical  problems which depend on fitting functions or in more specific fields like computer vision for problems related shape detection. 

For a brief description of the potential of this program and the difference between related strategies in the literature, like least square, consider the following dataset

``A=[-2.0  5.0; 
  -1.5  3.25;
  -1.0  2.0 ;
  -0.5  1.25;
   0.0  1.0 ;
   0.5  2.55;
   1.0  2.0 ;
   1.5  3.25;
   2.0  5.0 ;]``

Let us assume that we wish to find a fitting function given by parameters `x[1], x[2]` and the model

``f(t)=x[1]*t^2+x[2]``

where `f(t[i])≈A[i,2]` in the least square sense. When we apply a least square strategy we obtain

<img src="/docs/src/assets/figure_least_square_example.png" height="80%" width="80%"/>


and when we apply RAFF strategy we get

<img src="/docs/src/assets/figure_RAFF_example.png" height="80%" width="80%"/>

Note that RAFF detects (automatically) the model `f(t)=1.0*t^2 +1.0` that is, RAFF detects the point `(0.5,2.55)` like an outlier and excludes the influence of this point for determine `f(t)`. 
Consequently, RAFF detects the number of possible outliers and in the same time, which their are. This fact distinguishes RAFF from other techniques like [1] or classical least square.

Other useful examples are found in the documentation and also in directory [`examples/`](examples/).

## Installation
`RAFF` is a package registered in the `General repository`. Consequently, in order to install  `RAFF`, just type in `REPL`:

```julia
julia> ] add RAFF
```

## Advanced

Several files with test problems are provided inside directory
[`test/test_problems`](test/test_problems). We also provide some
scripts to create new random test problems, visualize and solve them
using `RAFF`. Those files can be found in
[`test/scripts`](test/scripts). For a detailed explanation of each
script file, please refer to the
[documentation](https://fsobral.github.io/RAFF.jl/stable/advanced/).

## Contributing

Any suggestions or contributions to `RAFF` are very welcome. `RAFF` is
a open-source software available through the MIT license, so users are
free to create their own versions of the software. Users are specially
encouraged to participate on the project in the following ways:

  - Modifications to the code can be submitter as pull requests to the
	branch `master`.
  - Bugs, suggestions and questions can be raised as issues in the
    project's Github page.
  - General doubts can also be sent directly to `fncsobral at uem dot
    br`.

## Performance

The following table outlines the results found by the current version
of `RAFF` in the test set.

| Name | Dim. | N Points | N Outl. | Found | Correct | Time (s) | Status | Solution |
| ---- | ---- | -------- | ------- | ----- | ------- | -------- | ------ | -------- |
|    C1 |   4 |     10 |   2 |   2 |   2 |     0.0006 | 1 | [ 1.999e+00,  1.142e-02, -4.023e+00, -1.053e+01] |
|    C2 |   4 |     10 |   2 |   1 |   1 |     0.0005 | 1 | [ 1.575e+00,  2.635e+00,  3.434e+01, -7.821e+00] |
|    C3 |   4 |     10 |   2 |   2 |   1 |     0.0005 | 1 | [ 4.221e-01, -1.596e+01,  2.388e+01,  2.146e+02] |
|    C4 |   4 |     10 |   1 |   1 |   1 |     0.0005 | 1 | [ 2.701e+00,  2.556e+00, -5.914e+01, -1.021e+02] |
|    C5 |   4 |    100 |   1 |   9 |   1 |     0.0210 | 1 | [ 2.114e+00,  1.160e+00, -6.590e+00, -1.971e+01] |
|   CI1 |   3 |     25 |   2 |   2 |   2 |     0.0020 | 1 | [ 9.741e-01,  9.769e-01,  5.020e+00] |
|    E1 |   3 |     20 |   1 |   1 |   0 |     0.0263 | 1 | [ 3.271e+03, -9.884e+01,  5.507e+01] |
|    E2 |   3 |     10 |   1 |   2 |   0 |     0.0110 | 1 | [ 2.137e+03,  2.760e+03,  2.369e+02] |
|    E3 |   3 |    100 |   1 |   1 |   1 |     0.2297 | 1 | [ 9.590e+02,  9.632e+03,  3.940e-01] |
|    E4 |   3 |    100 |  10 |   1 |   1 |     0.2308 | 1 | [ 1.139e+03,  9.370e+03,  4.206e-01] |
|  LOG1 |   4 |     10 |   1 |   0 |   0 |     0.0042 | 1 | [ 1.479e+03,  2.201e+03,  6.611e+01, -2.075e+02] |

## Credits

* Francisco Sobral (Leader)
* Emerson Vitor Castelani
* Ronaldo Lopes
* Wesley Shirabayashi

The authors of this package were sponsored by **Fundação Araucária**,
project number 002/17 - 47223. If you would like to cite this package,
please use

> Castelani, E. V., Lopes, R., Shirabayashi, W., & Sobral,
> F. N. C. (2019). RAFF.jl: Robust Algebraic Fitting Function in
> Julia. *Journal of Open Source Software*,
> 4(39), 1385. https://doi.org/10.21105/joss.01385

The following paper describes the theory and several comparison tests

> @article{Castelani2021,
>    author = {Castelani, Emerson V. and Lopes, Ronaldo and Shirabayashi, Wesley V. I. and Sobral, Francisco N. C.},
>    doi = {10.1007/s10898-020-00970-4},
>    journal = {Journal of Global Optimization},
>    title = {{A robust method based on LOVO functions for solving least squares problems}},
>    url = {http://link.springer.com/10.1007/s10898-020-00970-4},
>    year = {2021}
> }

## References

[1] Andreani, R., Martínez, J. M., Martínez, L., & Yano, F. S. (2009). Low order-value 
optimization and applications. *Journal of Global Optimization*, 43(1), 1-22.
