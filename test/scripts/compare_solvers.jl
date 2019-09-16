# This file runs different solvers and strategies for fitting data
# subject to noise and outliers. The codes were used from known
# sources and from different languages using Julia's support to Python
# and C.
#
# The following algorithms/strategies are tested
#
# - From Scipy: Least Squares, L1, Huber and Cauchy loss functions
# - From Theia: RANSAC implementations

# In order to run Scipy tests, one should install
# --------------------------
#
# Scipy 1.3.0 and Numpy 1.17.0 in the usual Python3 environment
#
#     pip3 install scipy numpy
#
# PyCall in Julia
#
#     ENV["PYTHON"] = "/usr/bin/python3"
#     ] add PyCall

# In order to run Theia test, we need to call C++ libraries, so it is
# necessary to install Theia.
#
# - Install Ceres
# - Add FindCeres.cmake and others
# - Install a lot of libraries

using DelimitedFiles
using PyCall
using PyPlot
using RAFF
using Printf
using Random

# Load libraries for drawing solutions
include("draw.jl")
include("gen_circle.jl")

"""

    compare_fitting(;model_str="linear", kwargs...)

Run several tests from Scipy.optimize library for fitting data.

"""
function compare_fitting(;model_str="linear", MAXMS=1, kwargs...)
    
    n, modl, = RAFF.model_list[model_str]

    optimize = pyimport("scipy.optimize")

    open("/tmp/output.txt") do fp
        
        global N = parse(Int, readline(fp))
        
        global data = readdlm(fp)
        
    end

    # SEEDMS is an argument for RAFF, the seed for multi-start strategy
    SEEDMS = 123456789
    
    # Create objective function for Python calls
    f = let

        x = copy(data[:, 1])

        y = copy(data[:, 2])

        c = zip(x, y)
    
        # z = zeros(length(x))
    
        θ -> map(t->modl(t[1], θ) .- t[2], c)
    end

    draw_problem(data)

    t = minimum(data[:, 1]):0.01:maximum(data[:, 1])

    # Run all fitting tests from Scipy
    
    for (i, (loss, line)) in enumerate(zip(
        ["linear", "soft_l1", "huber", "cauchy", "arctan"],
        ["--", "-.", ":", "-", "--"]
    ))
    
        seedMS = MersenneTwister(SEEDMS)

        fbest = Inf
        sbest = Nothing
        
        tm = @elapsed for i = 1:MAXMS

            initguess = randn(seedMS, Float64, n)

            sol = optimize.least_squares(f, initguess, loss=loss)

            if sol["cost"] < fbest

                fbest = sol["cost"]
                sbest = sol

            end
            
        end

        @printf("%10s & %8.4f & ", loss, tm)
        
        for j = 1:n

            @printf("\$%15.5f\$, ", sbest["x"][j])

        end

        @printf(" \\\\\n")

        modl2 = (x) -> modl(x, sbest["x"])

        PyPlot.plot(t, modl2.(t), color=PyPlot.cm."Set1"(i/9.0),
                    linestyle=line, label=loss)

    end

    # Run RAFF
    
    initguess = zeros(Float64, n)

    rsol, tm = @timed raff(modl, data[:, 1:end - 1], n; kwargs..., initguess=initguess,
                           SEEDMS=SEEDMS, MAXMS=MAXMS)

    @printf("%10s & %8.4f & ", "RAFF.jl", tm)

    for j = 1:n

        @printf("\$%15.5f\$, ", rsol.solution[j])

    end

    @printf(" \\\\\n")

    modl2 = (x) -> modl(x, rsol.solution)

    PyPlot.plot(t, modl2.(t), color=PyPlot.cm."Set1"(2.0/9.0),
                linestyle="-", label="RAFF")

    PyPlot.legend(loc="best")

    PyPlot.show()

    PyPlot.savefig("/tmp/scipy.png", DPI=150)
    
end

function compare_circle_fitting(;kwargs...)
    
    n, modl, = RAFF.model_list["circle"]

    optimize = pyimport("scipy.optimize")

    open("/tmp/output.txt") do fp
        
        global N = parse(Int, readline(fp))
        
        global data = readdlm(fp)
        
    end

    f = let

        l, = size(data)
        
        x = copy(data[:, 1:2])

        x = transpose(x)
        
        θ -> [modl(@view(x[:, i]), θ) for i = 1:l]
    end

    draw_circle_sol(data)

    t = 0:0.1:2.1 * π
    
    ptx = (α, ρ, d) -> ρ * cos(α) + d[1]
    pty = (α, ρ, d) -> ρ * sin(α) + d[2]    

    # Run all fitting tests from Scipy
    
    for (i, (loss, line)) in enumerate(zip(
        ["linear", "soft_l1", "huber", "cauchy", "arctan"],
        ["--", "-.", ":", "-", "--"]
    ))
    
        initguess = ones(Float64, n)

        sol, tm, = @timed optimize.least_squares(f, initguess, loss=loss)

        println("$(loss) -> $(tm)")

        fSol = sol["x"]
        
        modl1x = (α) -> ptx(α, fSol[3], fSol[1:2])
        modl1y = (α) -> pty(α, fSol[3], fSol[1:2])

        PyPlot.plot(modl1x.(t), modl1y.(t), color=PyPlot.cm."Set1"(i/9.0),
                    linestyle=line, label=loss)

    end

    # Run RAFF
    
    initguess = ones(Float64, n)

    rsol, tm = @timed raff(modl, data[:, 1:end - 1], n; kwargs..., initguess=initguess)

    println("RAFF -> $(tm)")
    println(rsol)
    
    fSol = rsol.solution
        
    modl1x = (α) -> ptx(α, fSol[3], fSol[1:2])
    modl1y = (α) -> pty(α, fSol[3], fSol[1:2])

    PyPlot.plot(modl1x.(t), modl1y.(t), color=PyPlot.cm."Set1"(2.0/9.0),
                linestyle="-", label="RAFF")

    PyPlot.legend(loc="best")

    PyPlot.show()

end
