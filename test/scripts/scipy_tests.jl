using DelimitedFiles
using PyCall
using PyPlot
using RAFF

include("draw.jl")

function compare_fitting(;model_str="linear", kwargs...)
    
    n, modl, = RAFF.model_list[model_str]

    optimize = pyimport("scipy.optimize")

    open("/tmp/output.txt") do fp
        
        global N = parse(Int, readline(fp))
        
        global data = readdlm(fp)
        
    end

    f = let

        x = copy(data[:, 1])

        y = copy(data[:, 2])

        c = zip(x, y)
    
        # z = zeros(length(x))
    
        θ -> map(t->modl(t[1], θ) .- t[2], c)
    end

    draw_problem(data)

    t = minimum(data[:, 1]):0.01:maximum(data[:, 1])

    for (i, (loss, line)) in enumerate(zip(
        ["linear", "soft_l1", "huber", "cauchy", "arctan"],
        ["--", "-.", ":", "-", "--"]
    ))
    
        initguess = zeros(Float64, n)

        sol, tm, = @timed optimize.least_squares(f, initguess, loss=loss)

        println("$(loss) -> $(tm)")

        modl2 = (x) -> modl(x, sol["x"])

        PyPlot.plot(t, modl2.(t), color=PyPlot.cm."Set1"(i/9.0),
                    linestyle=line, label=loss)

    end

    initguess = ones(Float64, n)

    rsol, tm = @timed raff(modl, data[:, 1:end - 1], n; kwargs..., initguess=initguess)

    println("RAFF -> $(tm)")
    println(rsol)

    modl2 = (x) -> modl(x, rsol.solution)

    PyPlot.plot(t, modl2.(t), color=PyPlot.cm."Set1"(2.0/9.0),
                linestyle="-", label="RAFF")

    PyPlot.legend(loc="best")

    PyPlot.show()

end
