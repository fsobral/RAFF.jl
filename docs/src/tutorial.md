# Tutorial

```@setup docrepl
```
## Installation

This package is supported just for Julia version 1.0. Consequently, 
it uses package 3.0. Currently `RAFF` is registered in [General Julia Registers](https://github.com/JuliaRegistries), so the 
package can be installed using the Julia package manager.
From the Julia REPL, type `]` to enter into Pkg REPL mode and run:

```
pkg> add RAFF
```

## Basic usage

Just to illustrate the potential and basic usage of `RAFF`, let us consider
the following data set given by an array ``A``:

```math
A=\left[ \begin{array}{cc}
 -2.0 &  5.0 \\
 -1.5 &  3.25\\
 -1.0 &  2.0 \\
 -0.5 &  1.25\\
  0.0 &  1.0 \\
  0.5 &  1.25\\
  1.0 &  2.0 \\
  1.5 &  3.25\\
  2.0 &  5.0 \\
\end{array}\right]
```

Let's suppose the first column of ``A`` as an experimental measure with 
result given by second column. It is easy to see that the fitting 
function in this case is accurate and given by 

```math
\phi(t)=t^2 +1
```

Now let's perturb one result of second column of ``A``. For example, 
consider ``A_{6,2} = 2.55``. Assuming the model for fitting given by
```math
\varphi(x,t)=x_1 t^2 +x_2 
```
we have by classical least squares as result `x = [0.904329,
1.3039]`. On the other hand, when we consider the `RAFF` algorithm we
obtain the correct answer `x=[1.0, 1.0]`. Moreover, we also have a
list of the possible outliers.

In order to run `RAFF` algorithm we need to setup 
```@repl docrepl
using RAFF
``` 
and define the data set and model:

```@repl docrepl
A=[-2.0  5.0; 
  -1.5  3.25;
  -1.0  2.0 ;
  -0.5  1.25;
   0.0  1.0 ;
   0.5  2.55;
   1.0  2.0 ;
   1.5  3.25;
   2.0  5.0 ;]

model(x, t) = x[1] * t[1]^2 + x[2]
```

After that, we can run method [`raff`](@ref):

```@repl docrepl
raff(model, A, 2)
```
The number `2` above is the number of variables in model, i. e., the
number of parameters to adjust in the model. The output is a
[`RAFFOutput`](@ref) type. For example to access only the parameters of
solution, we can use

```@repl docrepl
output = raff(model, A, 2)
output.solution
```

Note that `RAFF` algorithm detects and ignores possible outliers. In
order to see which points are outliers, we can access the `outliers`
attribute.

```@repl docrepl
output.outliers
```

More details about `RAFFOutput` type and other options can be obtained
in [API section](api.md).


By default `RAFF` uses automatic differentiation, more specifically
[ForwardDiff package](https://github.com/JuliaDiff/ForwardDiff.jl). But is possible
to call `RAFF` methods with gradient vector of model. For
example, considering the above example, we have,

```math
\nabla \varphi(x, t) = [t^2, 1].
```
Programming this gradient and run raff we have

```@repl docrepl
gmodel(x, t, g)=begin
   g[1] = t[1]^2
   g[2] = 1.0
   return g
end

raff(model, gmodel, A, 2)
```

Preliminary tests have shown that the use of explicit derivatives is
10 times faster than automatic differentiation.

## Multivariate models

`RAFF` supports the use of multivariate fitting functions to data sets
of different dimensions. To illustrate how this works, consider the
following example:

```@repl docrepl
data = [1.0 1.0   2.0
        0.0 0.0   4.0
        7.0 1.5  -4.5
        2.0 2.0 -17.0 # outlier
        0.0 8.6  -4.6]
```
and the following model

```@repl docrepl
model(x, t) = x[1] * t[1] + x[2] * t[2] + x[3]            
```

Note that this model has two variables ``(t_1, t_2)``. Naturally, this
problem has one outlier (`data[4,:]`), so there are 4 trust
points. Let's run `RAFF` and check the answer.

```@repl docrepl
output = raff(model, data, 3)
```

The right answer is `[- 1.0, - 1.0, 4.0]`. As we can note, `RAFF` get
a good fit for the data set. Handling the output follows the same
pattern as the one-dimensional case.

In order to get improvements in processing time, we can code the
gradient vector of model too:

```@repl docrepl
gmodel(x, t, g) = begin 
    g[1] = t[1]
    g[2] = t[2]
    g[3] = 1.0
    return g
end
```
```@repl docrepl
output = raff(model,data, 3)
```


## Changing some options

Naturally, `RAFF` has options like precision of gradient stopping criteria and initial guess. 

```@repl docrepl
output = raff(model,data, 3;initguess=[0.5,0.5,0.5],ε=1.0e-4)
```

`RAFF` is based on an optimization method. In this way, it is subject to
stopping at stationary points that are not global minimizers. For this
reason, heuristics were implemented to find global minimizers. Such
heuristics depend on random number generation. So, if you want to run
tests with more reliability this can be a useful strategy. To define
in `RAFF`, say, 1000 different starting points, is enough to redefine
the keyword argument `MAXMS`.

```@repl docrepl
output = raff(model, data, 3; MAXMS=1, initguess=[0.5,0.5,0.5], ε=1.0e-10)
```

In the above example, we have also changed the starting point for the
method. Also, the stopping criterion was changed to ``10^{-10}``,
which means high accuracy when solving the subproblems. See [`RAFF`
API](api.md#RAFF) for all the possible options that can be used.


## Parallel running

`RAFF` can be run in a parallel or distributed environment, using the
[Distributed](https://docs.julialang.org/en/v1.0/stdlib/Distributed/)
package and function [`praff`](@ref). Let's use `praff` to solve the
same problem from the beginning. First, the Distributed package has to
be loaded and the number of workers has to be added. It is also
possible to add the address of other machines.

```
using Distributed

addprocs(3) # Add 3 worker processes
```

This step can be replaced if Julia is initialized with the `-p`
option

```
julia -p 3
```

Now we have to load [`RAFF`](@ref Overview) and the fit function in all
workers:

```
@everywhere using RAFF

@everywhere function gmodel!(x, t, g)
    g[1] = t[1]^2
    g[2] = 1.0
end

@everywhere function model(x, t)
   x[1] * t[1]^2 + x[2]
end
```

then, we call [`praff`](@ref) to solve the problem (note that we do
not need to send the `A` matrix to all workers, since it will be
automatically sent by `praff`).

```
A=[-2.0  5.0;
  -1.5  3.25;
  -1.0  2.0 ;
  -0.5  1.25;
   0.0  1.0 ;
   0.5  2.55;
   1.0  2.0 ;
   1.5  3.25;
   2.0  5.0 ;];

n = 2

output = praff(model, gmodel!, A, n)
RAFFOutput(1, [1.0, 0.999996], 6, 8, 4.0205772365906425e-11, [6])
```

The true effectiveness of parallelism occurs when option `MAXMS` is
set, which changes the number of random initial points that are tried
for each subproblem solved. Better solutions can be achieved with
higher values of `MAXMS`

```
n = 2

output = praff(model, gmodel!, A, n; MAXMS=1000)
RAFFOutput(1, [1.0, 1.0], 7, 8, 5.134133698545651e-13, [6])
```