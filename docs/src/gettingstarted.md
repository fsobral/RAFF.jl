## Getting Started

## Installation

This package is supported just for Julia version 1.0. Consequently, 
it uses package 3.0. Currently RAFF ins't in 
[Metadata.jl](https://github.com/JuliaLang/METADATA.jl), so the 
package can be installed with the Julia package manager.
From the Julia REPL, type `]` to enter into Pkg REPL mode and run:

```
pkg> add https://github.com/fsobral/RAFF.jl#dev
```

## Basic usage

Just to ilustrate the potential and basic usage of RAFF, let us consider
the following dataset given by an array ``A``:

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

Lets suppose the first column of ``A`` as an experimental measure with 
result given by second column. It is easy to see that the fitting 
function in this case is accurate and given by 

```math
\phi(t)=t^2 +1
```

Now lets perturb one result of second column of ``A``. For example, 
consider ``A_{6,2}=2.55``. Assuming the model for fitting given by

```math
\varphi(x,t)=x_1 t^2 +x_2 
```
we have by classical least square as result `x=[0.04333,2.8111]`. But when we consider RAFF algorithm we obtain the correct answer `x=[1.0,1.0]`. Moreover, we have the number of possible outliers and which they are.

In order to get run RAFF algorithm we need to setup 
```@repl 1
using RAFF
``` 
and define the dataset and model:

```@repl 1
A=[-2.0  5.0; 
  -1.5  3.25;
  -1.0  2.0 ;
  -0.5  1.25;
   0.0  1.0 ;
   0.5  1.25;
   1.0  2.0 ;
   1.5  3.25;
   2.0  5.0 ;]

model(x,t)=x[1]*t^2+x[2]
```

After that, we can run raff:

```@repl 1
raff(model,A,2)
```



