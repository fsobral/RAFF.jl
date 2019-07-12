export RAFFOutput

"""

This type defines the output file for the RAFF algorithm.

    RAFFOutput(status::Int, solution::Vector{Float64}, iter::Int,
               p::Int, f::Float64, outliers::Vector{Int})

where

  - `status`: is 1 if converged and 0 if not
  - `solution`: vector with the parameters of the model
  - `iter`: number of iterations up to convergence
  - `p`: number of trusted points
  - `f`: the residual value
  - `outliers`: the possible outliers detected by the method, for the
    given `p`


    RAFFOutput()

Creates a null version of output, equivalent to `RAFFOutput(0, [], -1, 0, Inf, [])`

    RAFFOuput(p::Int)

Creates a null version of output for the given `p`.

"""
struct RAFFOutput

    status :: Int
    solution :: Vector{Float64}
    iter :: Int
    p :: Int
    f :: Float64
    outliers :: Vector{Int}
    
end

RAFFOutput(p) = RAFFOutput(0, [], -1, p, Inf, [])

RAFFOutput() = RAFFOutput(0)

# Overload the == for RAFFOutput

import Base.==

==(a::RAFFOutput, b::RAFFOutput) = ((a.status == b.status) &&
     (a.solution == b.solution) && (a.iter == b.iter) &&
     (a.p == b.p) && (a.f == b.f) && (a.outliers == b.outliers))

import Base.show

function show(io::IO, ro::RAFFOutput)
	print(io,"** RAFFOutput ** \n")
	print(io,"Status (.status) = $(ro.status) \n")
	print(io,"Solution (.solution) = $(ro.solution) \n")
	print(io,"Number of iterations (.iter) = $(ro.iter) \n")
	print(io,"Number of trust points (.p) = $(ro.p) \n")
	print(io,"Objective function value (.f) = $(ro.f) \n")
	print(io,"Index of outliers (.outliers)= $(ro.outliers)\n")
end