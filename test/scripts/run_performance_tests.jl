# This script runs all the tests and create a table for comparison.

using RAFF
using Printf
using DelimitedFiles

fpath = "../test_problems/"

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
    
    @printf(tfp, """
| Name | Dim. | N Points | N Outl. | Found | Correct | Time (s) | Status | Solution |
| ---- | ---- | -------- | ------- | ----- | ------- | -------- | ------ | -------- |
""")
        
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

        # Save problem data

        name = match(r"[^.]*", fname).match

        @printf(tfp, "| %5s | %3d | %6d | %3d | ", name, n, np, length(outliers))

        nexact = 0

        for i in rsol.outliers

            (i in outliers) && (nexact += 1)

        end

        @printf(tfp, "%3d | %3d | %10.4f | %d | ",
                length(rsol.outliers), nexact, t, rsol.status)

        # Save best solution
        
        @printf(tfp, "[")

        for i = 1:n - 1

            @printf(tfp, "%10.3e, ", rsol.solution[i])

        end

        @printf(tfp, "%10.3e] |\n", rsol.solution[n])

    end

    close(tfp)

end
