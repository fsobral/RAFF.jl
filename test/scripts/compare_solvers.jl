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

function run_comparative_fitting()

    large_number = 179424673

    for (model_str, sol) in [ ("linear", [-200.0, 1000.0]), ("cubic", [0.5, -20.0, 300.0, 1000.0]),
                              ("expon", [5000.0, 4000.0, 0.2]),
                              ("logistic", [6000.0, -5000, -0.2, -3.7]) ]

        for (np, p) in [(10, 9), (10, 8), (100, 99), (100, 90)]

            n, model, = RAFF.model_list[model_str]

            # Define seed for this run. The same seed for all instances.
            Random.seed!(large_number + 300)

            data, = generate_noisy_data(model, n, np, p, sol, (1.0, 30.0), std=200.0)

            cla()
            
            compare_fitting(1, data; model_str=model_str, MAXMS=200,
                            outimage="/tmp/$(model_str)_$(n)_$(np).eps")

        end

    end
    
end

"""

    ls_measure(modl, n, data)

Compute the "Least Square" measure, used to compare different
algorithm.

"""

function ls_measure(modl, n, data)

    sr = 0.0

    np, = size(data)
        
    for i = 1:np

        (data[i, n + 2] == 0.0) && (sr += (modl(@view(data[i, 1:n])) - data[i, n + 1])^2)

    end

    return sr

end

"""

    compare_fitting(;model_str="linear", kwargs...)

Run several tests from Scipy.optimize library for fitting data.

"""
function compare_fitting(N=Nothing, data=Nothing;model_str="linear", MAXMS=1, outimage="/tmp/scipy.png",
                         kwargs...)
    
    n, modl, = RAFF.model_list[model_str]

    optimize = pyimport("scipy.optimize")

    if (N == Nothing) || (data == Nothing)

        open("/tmp/output.txt") do fp
        
            N = parse(Int, readline(fp))
        
            data = readdlm(fp)
        
        end

    end

    # SEEDMS is an argument for RAFF, the seed for multi-start strategy
    SEEDMS = 123456789

    # PyPlot configuration
    PyPlot.rc("font", family="Helvetica")

    PyPlot.rc("font", size=10)
    
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

    tmpv = @view(data[:, 3])
    
    @printf("\\multirow{5}{*}{\$(%10s, %4d, %4d)\$} \n", model_str, length(tmpv), count(tmpv .== 0.0))
    
    # Run all fitting tests from Scipy. Not using 'arctan', due to its
    # terrible results.
    
    for (i, (loss, line)) in enumerate(zip(
        ["linear", "soft_l1", "huber", "cauchy"],
        ["--", "-.", ":", "-"]
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

        # This is a tentative of comparing different algorithms.  We
        # compute the sum of the squared differences only for
        # NON-outliers.
        
        modl2 = (x) -> modl(x, sbest["x"])

        fmeas = ls_measure(modl2, N, data)
                           
        @printf("  & %10s & %10.3e & %8.4f & %6d & ", loss, fmeas, tm, sbest["nfev"])
        
        for j = 1:n

            @printf("\$%15.5f\$, ", sbest["x"][j])

        end

        @printf(" \\\\\n")

        PyPlot.plot(t, modl2.(t), color=PyPlot.cm."Set1"(i/9.0),
                    linestyle=line, label=loss)

    end

    # Run RAFF
    
    initguess = zeros(Float64, n)

    rsol, tm = @timed raff(modl, data[:, 1:end - 1], n; kwargs..., initguess=initguess,
                           SEEDMS=SEEDMS, MAXMS=MAXMS)

    modl2 = (x) -> modl(x, rsol.solution)

    fmeas = ls_measure(modl2, N, data)
                           
    @printf("  & %10s & %10.3e & %8.4f & %6d & ", "RAFF.jl", fmeas, tm, 0)
        
    for j = 1:n

        @printf("\$%15.5f\$, ", rsol.solution[j])

    end

    @printf(" \\\\\n")

    PyPlot.plot(t, modl2.(t), color=PyPlot.cm."Set1"(2.0/9.0),
                linestyle="-", label="RAFF")

    PyPlot.legend(loc="best")

    PyPlot.show()

    PyPlot.savefig("/tmp/comp_scipy.eps", DPI=600, bbox_inches="tight")

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
