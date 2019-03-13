# This script runs all the tests and create a table for comparison.

using RAFF
using Printf
using DelimitedFiles

fpath = "../test/test-problems/"

getType = begin

    d = Dict(
        "C" => "cubic",
        "L" => "linear",
        "LOG" => "logistic",
        "CI" => "circle",
        "E" => "expon"
    )

    (s) -> d[s]

end

function run_tests(;kwargs...)

    tfp = open("/tmp/table.csv", "w")
    
    for fname in readdir(fpath)

        occursin("sol", fname) && continue

        m = match(r"[^1-9]*", fname)

        (m == nothing) && continue

        mtype = getType(m.match)

        open(fpath * fname) do fp

            global N = parse(Int, readline(fp))

            global data = readdlm(fp)

        end

        n, model, = RAFF.model_list[mtype]

        np, = size(data)
        
        # Retrieve the true outliers
        outliers = findall(data[:, end] .== 1.0)

        # Run RAFF
        rsol, t, = @timed raff(model, data[:, 1:end - 1], n; kwargs...)

        name = match(r"[^.]*", fname).match
        
        @printf(tfp, "%5s\t%6d\t%3d\t[", name, np, length(outliers))

        for i = 1:n - 1

            @printf(tfp, "%10.3e, ", rsol.solution[i])

        end

        @printf(tfp, "%10.3e]\t", rsol.solution[n])

        nexact = 0

        for i in rsol.outliers

            (i in outliers) && (nexact += 1)

        end

        @printf(tfp, "%3d\t%3d\t%10.4f\t%d\n",
                length(rsol.outliers), nexact, t, rsol.status)

    end

    close(tfp)

end
