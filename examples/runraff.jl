using RAFF
using DelimitedFiles
using Printf

# This file loads and runs all test files. It prints the obtained and
# expected outputs.

# Iterate over a list of problems and solutions
for prob in eachline("list.dat")

    # Ignore blank lines
    (length(strip(prob)) == 0) && continue
    
    dname, sname = split(prob)

    # Data file
    data = readdlm(dname)[:, [1, 2]]

    # Solution file
    fsol = open(sname, "r")

    # Number of parameters
    n = Meta.parse(readline(fsol))

    # Solution vector
    tsol = eval(Meta.parse(readline(fsol)))

    # Model function to fit data
    model = eval(Meta.parse(readline(fsol)))

    close(fsol)

    # Call RAFF
    s, x, iter, p = raff(model, data, n)

    # Output information
    @printf("\nFile: %s\nSolved:\n\tStats: %d\n\tIters: %d\n\tp: %d\n",
            dname, s, iter, p)
    println("Found:   ", x)
    println("Expected:", tsol)
    
end 
