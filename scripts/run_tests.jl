using RAFF
using Printf
using DelimitedFiles

fpath = "../test/test-problems/"
fname = "CI1"
mtype = "circle"

open(fpath * fname * ".txt") do fp
        
    global N = parse(Int, readline(fp))
        
    global data = readdlm(fp)
        
end

n, model, = RAFF.model_list[mtype]

# Retrieve the true outliers
outliers = findall(data[:, end] .== 1.0)

# Run RAFF
rsol, t, = @timed raff(model, data[:, 1:end - 1], n)

open("/tmp/table.csv", "a") do fp

    @printf(fp, "%5s\t[", fname)

    for i = 1:n - 1

        @printf(fp, "%10.3e, ", rsol.solution[i])

    end

    @printf(fp, "%10.3e]\t", rsol.solution[n])

    nexact = 0

    for i in rsol.outliers

        (i in outliers) && (nexact += 1)

    end
    
    @printf(fp, "%3d\t%3d\t%10.4f\n", length(rsol.outliers), nexact, t)

end
