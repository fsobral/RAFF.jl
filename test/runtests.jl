using RAFF

using Test
using DelimitedFiles
using Distributed
using Random
using SharedArrays
using Logging
using Base.CoreLogging

global_logger(ConsoleLogger(stdout, Logging.Debug))

include("test_parallel.jl")

# include("test_utils.jl")

# include("test_raff.jl")

# include("test_integration.jl")
