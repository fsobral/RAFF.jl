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

# In order to run Theia test, we have developed a C++ program that
# runs RANSAC together with CERES solver from Google. The code and
# more details are available in `../theia` directory.

using DelimitedFiles
using PyCall
using PyPlot
using RAFF
using Printf
using Random
using Base.GC

# Load libraries for drawing solutions
include("draw.jl")
include("gen_circle.jl")

"""

    run_comparative_fitting()

Run all the regular tests with one-dimensional problems and all the
solvers.

"""
function run_comparative_fitting()

    large_number = 179424673

    # First call. Just compilation.
    compare_fitting(1, [1.0 1.0 0.0; 2.0 2.0 0.0])

    for (model_str, sol) in [ ("linear", [-200.0, 1000.0]), ("cubic", [0.5, -20.0, 300.0, 1000.0]),
                              ("expon", [5000.0, 4000.0, 0.2]),
                              ("logistic", [6000.0, -5000, -0.2, -3.7]) ]

        for (np, p) in [(10, 9), (10, 8), (100, 99), (100, 90)]

            n, model, mstr = RAFF.model_list[model_str]

            # Define seed for this run. The same seed for all instances.
            Random.seed!(large_number + 300)

            generate_test_problems("/tmp/output.txt", "/tmp/sol.txt", model, mstr, n, np, p;
                                   θSol=sol, x_interval=(1.0, 30.0), std=100.0)

            cla()
            
            compare_fitting(model_str=model_str, MAXMS=100, theia_ms=10,
                            outimage="/tmp/$(model_str)_$(np)_$(p).png")

        end

    end

    # Run cluster tests

    for (model_str, sol) in [ ("linear", [-200.0, 1000.0]), ("cubic", [0.5, -20.0, 300.0, 1000.0]),
                              ("expon", [5000.0, 4000.0, 0.2]),
                              ("logistic", [6000.0, -5000, -0.2, -3.7]) ]

        for (np, p) in [(10, 8), (100, 90)]

            n, model, mstr = RAFF.model_list[model_str]

            # Define seed for this run. The same seed for all instances.
            Random.seed!(large_number + 300)

            generate_test_problems("/tmp/output.txt", "/tmp/sol.txt", model, mstr, n, np, p,
                                   (1.0, 30.0), (5.0, 10.0); θSol=sol, std=100.0)

            cla()
            
            compare_fitting(model_str=model_str, MAXMS=100, theia_ms=10,
                            outimage="/tmp/$(model_str)_$(np)_$(p)_cluster.png")

        end

    end
    
end

"""

    ls_measure(modl, n, data)

Compute the "Least Square" measure, used to compare different
algorithms. By definition is the square root of the sum of the squared
difference between the value predicted by the model and the actual
value only for *non-outliers*.

"""

function ls_measure(modl, n, data)

    sr = 0.0

    np, = size(data)
        
    for i = 1:np

        (data[i, n + 2] == 0.0) && (sr += (modl(@view(data[i, 1:n])) - data[i, n + 1])^2)

    end

    return sqrt(sr)

end

"""

    compare_fitting(;model_str="linear", kwargs...)

Run several tests from Scipy.optimize library for fitting data.

"""
function compare_fitting(N, data; kwargs...)

    open("/tmp/output.txt", "w") do fp

        @printf(fp, "%d\n", N)

        l, c = size(data)
        
        for i = 1:l

            for j = 1:c

                @printf(fp, "%f ", data[i, j])

            end

            @printf(fp, "\n")

        end
        
    end

    compare_fitting(kwargs...)

end


function compare_fitting(filename="/tmp/output.txt"; model_str="linear", MAXMS=1, outimage="/tmp/scipy.png",
                         theia_str_args="2000 -ft 0.1", theia_ms=10, kwargs...)

    n, modl, = RAFF.model_list[model_str]

    optimize = pyimport("scipy.optimize")

    # SEEDMS is an argument for RAFF, the seed for multi-start strategy
    SEEDMS = 123456789

    # PyPlot configuration
    PyPlot.rc("font", family="Helvetica")

    PyPlot.rc("font", size=10)

    # Open file
    data = Array{Float64, 2}(undef, 0, 0)
    
    N = 1
    
    open(filename) do fp
        
        N = parse(Int, readline(fp))
        
        data = readdlm(fp)
        
    end

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

    open("table.txt", "a") do fp

        @printf(fp, "\\multirow{5}{*}{\$(%10s, %4d, %4d)\$} \n",
                model_str, length(tmpv), count(tmpv .== 0.0))

    end

    # Run all fitting tests from Scipy. Not using 'arctan', due to its
    # terrible results.

    for (i, (loss, line)) in enumerate(zip(
        ["linear", "soft_l1", "huber", "cauchy"],
        ["--", "-.", ":", "-"]
    ))
    
        seedMS = MersenneTwister(SEEDMS)

        fbest = Inf
        sbest = Nothing
        nfev  = 0

        # Eliminate effects of garbage collector
        GC.gc()
        
        tm = @elapsed for i = 1:MAXMS

            initguess = randn(seedMS, Float64, n)

            sol = optimize.least_squares(f, initguess, loss=loss)

            nfev += sol["nfev"]

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

        open("table.txt", "a") do fp

            @printf(fp, "  & %10s & %10.3e & %8.4f & %8d & ", loss, fmeas, tm, nfev)

            for j = 1:n

                @printf(fp, "\$%15.5f\$, ", sbest["x"][j])

            end

            @printf(fp, " \\\\\n")

        end

        PyPlot.plot(t, modl2.(t), color=PyPlot.cm."Set1"(i/9.0),
                    linestyle=line, label=loss)

    end

    # Run Theia

    fname = "/tmp/" * Random.randstring(10) * ".txt"

    run(pipeline(`../theia/ransac $(model_str) $(theia_str_args) -ms $(theia_ms) -f $filename`, fname))

    theia = DelimitedFiles.readdlm(fname)

    tt = minimum(data[:, 1]):1:maximum(data[:, 1])

    tm, = size(theia)

    for i = 1:Int(tm / theia_ms)

        b, bpos = findmin(@view(theia[(i - 1) * theia_ms + 1:i * theia_ms, 4]))

        total_time = sum(@view(theia[(i - 1) * theia_ms + 1:i * theia_ms, 2]))

        tvec = @view theia[(i - 1) * theia_ms + bpos, :]
        
        modl2 = (x) -> modl(x, tvec[end - n + 1:end])

        fmeas = ls_measure(modl2, N, data)

        open("table.txt", "a") do fp

            @printf(fp, "  & %10s & %10.3e & %8.4f & %8d & ", tvec[1], fmeas, total_time, 0)

            for j = 1:n

                @printf(fp, "\$%15.5f\$, ", tvec[end - n + j])

            end

            @printf(fp, " \\\\\n")

        end

        PyPlot.plot(tt, modl2.(tt), color=PyPlot.cm."Set1"((5 + i)/9.0),
                    linestyle="-", marker=3 + i, label=tvec[1])

    end

    # Run RAFF
    
    initguess = zeros(Float64, n)

    # Eliminate effects of garbage collector
    GC.gc()
        
    rsol, tm = @timed praff(modl, data[:, 1:end - 1], n; kwargs..., initguess=initguess,
                            SEEDMS=SEEDMS, MAXMS=MAXMS)

    modl2 = (x) -> modl(x, rsol.solution)

    fmeas = ls_measure(modl2, N, data)

    open("table.txt", "a") do fp

        @printf(fp, "  & %10s & %10.3e & %8.4f & %8d & ", "RAFF.jl", fmeas, tm, rsol.nf)

        for j = 1:n

            @printf(fp, "\$%15.5f\$, ", rsol.solution[j])

        end

        @printf(fp, " \\\\\n")

    end

    PyPlot.plot(t, modl2.(t), color=PyPlot.cm."Set1"(2.0/9.0),
                linestyle="-", label="RAFF")

    PyPlot.legend(loc="best")

    PyPlot.show()

    PyPlot.savefig(outimage, DPI=600, bbox_inches="tight")

    PyPlot.savefig("/tmp/scipy.png", DPI=150)
    
end

function compare_circle_fitting(filename="/tmp/output.txt";MAXMS=1, outimage="/tmp/scipy.png",
                                theia_str_args="10 -ft 0.1", kwargs...)
    
    n, modl, = RAFF.model_list["circle"]

    optimize = pyimport("scipy.optimize")

    # SEEDMS is an argument for RAFF, the seed for multi-start strategy
    SEEDMS = 123456789

    # PyPlot configuration
    PyPlot.rc("font", family="Helvetica")

    PyPlot.rc("font", size=10)

    # Open file
    data = Array{Float64, 2}(undef, 0, 0)

    N = 1

    open(filename) do fp

        N = parse(Int, readline(fp))

        data = readdlm(fp)

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

    tmpv = @view(data[:, 4])

    @printf("\\multirow{5}{*}{\$(%10s, %4d, %4d)\$} \n", "circle", length(tmpv), count(tmpv .== 0.0))

    # Run all fitting tests from Scipy. Not using 'arctan', due to its
    # terrible results.
    
    for (i, (loss, line)) in enumerate(zip(
        ["linear", "soft_l1", "huber", "cauchy"],
        ["--", "-.", ":", "-"]
    ))
    
        initguess = ones(Float64, n)

        seedMS = MersenneTwister(SEEDMS)

        fbest = Inf
        sbest = Nothing
        nfev  = 0
        
        tm = @elapsed for i = 1:MAXMS

            initguess = 1.0 .+ randn(seedMS, Float64, n)

            sol = optimize.least_squares(f, initguess, loss=loss)

            nfev += sol["nfev"]

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

        @printf("  & %10s & %10.3e & %8.4f & %8d & ", loss, fmeas, tm, nfev)

        for j = 1:n

            @printf("\$%15.5f\$, ", sbest["x"][j])

        end

        @printf(" \\\\\n")

        modl1x = (α) -> ptx(α, sbest["x"][3], sbest["x"][1:2])
        modl1y = (α) -> pty(α, sbest["x"][3], sbest["x"][1:2])

        PyPlot.plot(modl1x.(t), modl1y.(t), color=PyPlot.cm."Set1"(i/9.0),
                    linestyle=line, label=loss)

    end

    # Run Theia

    fname = "/tmp/" * Random.randstring(10) * ".txt"

    multistart = MAXMS
    
    run(pipeline(`../theia/ransac circle $(theia_str_args) -ms $(MAXMS)`, fname))

    theia = DelimitedFiles.readdlm(fname)

    tt = 0:0.2:2.1 * π

    tm, = size(theia)

    for i = 1:Int(tm / multistart)

        b, bpos = findmin(@view(theia[(i - 1) * multistart + 1:i * multistart, 4]))

        total_time = sum(@view(theia[(i - 1) * multistart + 1:i * multistart, 2]))

        tvec = @view theia[(i - 1) * multistart + bpos, :]
        
        modl2 = (x) -> modl(x, tvec[end - n + 1:end])

        fmeas = ls_measure(modl2, N, data)

        @printf("  & %10s & %10.3e & %8.4f & %8d & ", tvec[1], fmeas, total_time, 0)

        for j = 1:n

            @printf("\$%15.5f\$, ", tvec[end - n + j])

        end

        @printf(" \\\\\n")

        modl1x = (α) -> ptx(α, tvec[end], tvec[end - 2:end - 1])
        modl1y = (α) -> pty(α, tvec[end], tvec[end - 2:end - 1])

        PyPlot.plot(modl1x.(tt), modl1y.(tt), color=PyPlot.cm."Set1"((5 + i)/9.0),
                    linestyle="-", marker=i, label=tvec[1])

    end

    # Run RAFF
    
    initguess = ones(Float64, n)

    rsol, tm = @timed praff(modl, data[:, 1:end - 1], n; kwargs..., initguess=initguess)

    fSol = rsol.solution
        
    modl2 = (x) -> modl(x, rsol.solution)

    fmeas = ls_measure(modl2, N, data)
                           
    @printf("  & %10s & %10.3e & %8.4f & %8d & ", "RAFF.jl", fmeas, tm, rsol.nf)
        
    for j = 1:n

        @printf("\$%15.5f\$, ", fSol[j])

    end

    @printf(" \\\\\n")
    
    modl1x = (α) -> ptx(α, fSol[3], fSol[1:2])
    modl1y = (α) -> pty(α, fSol[3], fSol[1:2])

    PyPlot.plot(modl1x.(t), modl1y.(t), color=PyPlot.cm."Set1"(2.0/9.0),
                linestyle="-", label="RAFF")

    PyPlot.legend(bbox_to_anchor=(0., 1.01, 1., .111), loc="lower left",
           ncol=3, mode="expand", borderaxespad=0., numpoints=1)

    PyPlot.show()

    PyPlot.savefig(outimage, DPI=600, bbox_inches="tight")

    PyPlot.savefig("/tmp/scipy.png", DPI=150)
    
end
