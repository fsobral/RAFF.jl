RAFF - Robust Algebraic Fitting Function
========================================

[![Build Status](https://travis-ci.com/fsobral/RAFF.jl.svg?token=AFD4xpyZzz3My5T6Sx1k&branch=master)](https://travis-ci.com/fsobral/RAFF.jl)
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://fsobral.github.io/RAFF.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://fsobral.github.io/RAFF.jl/dev)

This package implements a robust method to fit a given function
(described by its parameters) to input data. The method is based on
the LOVO algorithm [1] and also in a suitable voting strategy in order
automatically eliminate outliers.

The potential use of this package is in statistical  problems which depend on fitting functions or in more specific fields like computer vision for problems related shape detection. 

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

Note that RAFF detects (automatically) the model `f(t)=1.0*t^2 +1.0` that is, RAFF detects the point `(0.5,2.55)` like an outilier and excludes the influence of this point for determine `f(t)`. 
Consequently, RAFF detects the number of possible outliers and in the same time, which their are. This fact distinguishes RAFF from other techniques like [1] or classical least square.

## Credits

* Francisco Sobral (Leader)
* Emerson Vitor Castelani
* Ronaldo Lopes
* Wesley Shirabayashi

## References

[1] Andreani, R., Martínez, J. M., Martínez, L., & Yano, F. S. (2009). Low order-value 
optimization and applications. *Journal of Global Optimization*, 43(1), 1-22.