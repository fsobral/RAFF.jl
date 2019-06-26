using DelimitedFiles
using PyPlot
using RAFF

"""

    draw_problem(;raff_output=nothing, model_str="logistic",
                  datafile="/tmp/output.txt")

Draw the problem data. By default it is assumed that the model is
given by the logistic model and the data file is given in file
`/tmp/output.txt`.

If a [RAFFOutput](@ref) object is provided, then it plots the model
found and also the outliers (true and false positives).

Optional arguments:

  - `raff_output`: [RAFFOutput](@ref) object with the solution
    obtained
  - `model_str`: a string with the name of the model to be used to
    plot the solution. See [model_list](@ref) for details.
  - `datafile`: file containing the data.

"""
function draw_problem(;raff_output=nothing, model_str="logistic", datafile="/tmp/output.txt")

    fp = open(datafile, "r")

    N = parse(Int, readline(fp))

    M = readdlm(fp)

    close(fp)

    x = M[:, 1]
    y = M[:, 2]
    co = M[:, 3]

    true_outliers = findall(co .== 1)

    PyPlot.scatter(x[co .== 0.0], y[co .== 0.0], color=PyPlot.cm."Pastel1"(2.0/9.0),
                   marker="o", s=50.0, linewidths=0.2)

    PyPlot.scatter(x[co .== 1.0], y[co .== 1.0], color=PyPlot.cm."Pastel1"(2.0/9.0),
                   marker="^", s=50.0, linewidths=0.2, label="Outliers")

    if raff_output != nothing
        
        n, model, modelstr = RAFF.model_list[model_str]

        modl1 = (x) -> model(x, raff_output.solution)

        t = minimum(x):0.01:maximum(x)
        PyPlot.plot(t, modl1.(t), color=PyPlot.cm."Set1"(2.0/9.0))

        # Draw outliers found by RAFF
        
        true_positives = intersect(true_outliers, raff_output.outliers)
        false_positives = setdiff(raff_output.outliers, true_positives)

        PyPlot.scatter(x[false_positives], y[false_positives],
                       color=PyPlot.cm."Pastel1"(0.0/9.0), marker="o",
                       linewidths=0.2, edgecolors="k", s=50.0, label="False positives")
        
        PyPlot.scatter(x[true_positives], y[true_positives],
                       color=PyPlot.cm."Pastel1"(0.0/9.0), marker="^",
                       s=50.0, linewidths=0.2, edgecolors="k", label="Identified outliers")

    end

    # θSol = [1000.0, 5000.0, 0.2, 3.0]
    # PyPlot.plot(t, modl2.(t), "b--", label="Verdadeiro")

    PyPlot.legend(loc="best")

    PyPlot.show()
    
    PyPlot.savefig("/tmp/figure.png", DPI=100)

end
