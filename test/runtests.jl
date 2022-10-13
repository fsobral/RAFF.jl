using RAFF

using Test
using DelimitedFiles
using Distributed
using Random
using SharedArrays
using Logging
using LinearAlgebra

include("test_raff.jl")

include("test_solvers.jl")

include("test_smartqr.jl")

include("test_utils.jl")

include("test_generator.jl")

include("test_parallel.jl")

include("test_multivariate.jl")

include("test_integration.jl")
