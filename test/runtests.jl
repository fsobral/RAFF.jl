using RAFF

using Test
using DelimitedFiles
using Base.CoreLogging

# Disabling info and debug from logging
disable_logging(CoreLogging.Info)

include("test_utils.jl")

include("test_raff.jl")

include("test_integration.jl")
