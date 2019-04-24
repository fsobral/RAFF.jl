# Advanced usage

## Test Problems

In order to develop a robust code, test problems regarding data
fitting with outliers have been created. All the problems are
available in directory `test/test_problems`. All the problems have the
form

```math
\begin{array}{ccccc}
	k
	t_{11} & t_{12} & ... & t_{1k} & f_1\\
	... & & & & \\
	t_{np,1} & t_{np,2} & ... & t_{np,k} & f_{np}
\end{array}
```

where ``k`` is the number of *variables* of the model (not the number
of parameters!), ``np`` is the number of data points to be adjusted,
``t_{ij}`` are the values selected for parameter ``j`` in experiment
``i`` and ``f_i`` is the result obtained for experiment ``i``.

## Script files

During the development and testing of `RAFF` several scripts and
pieces of Julia code have been created. Those files are mostly related
to the automated generation of test problems and visualization of the
solutions. All those files are located in the `test/scripts`
directory.

We explain each file below, so maybe more advanced users can modify
and re-use the code to generate their own problems.

  - `calc_ratio.jl`: this script contains a function that generates
    several tests for the same model and solve each one with
    `RAFF`. Then it prints the ratio of outliers that have
    successfully been detected but the method.
  - `draw.jl`: This script draws the solution of a given problem and
    also its data, which was taken from a file. Examples of such files
    are located in `test/test_problems`.
  - `run_raff.jl`: this script simply loads a problem data from a
    file, selects a given model and runs `RAFF`. Examples of such
    files are located in `test/test_problems`.
  - `run_lmlovo.jl`: the same as before, but only for [`lmlovo`](@ref)
    function. Used mostly for testing.
  - `generate_fit_tests.jl`: script for generating random test problem
    files, using the pre-defined models given by
    [`RAFF.model_list`](@ref). This function cannot be called inside
    Julia, since it uses `ArgParse` package.
  - `gen_circle.jl`: specific script for generating random test problems
    related to the detection of circles in the plane. It also provides
    functions to draw the problem and the solution, which differ from
    the `draw.jl` script above.
  - `run_performance_tests.jl`: script for generating some performance
    tests, so we can compare different versions of RAFF.
