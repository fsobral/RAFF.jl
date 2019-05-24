using DelimitedFiles
using PyPlot
using RAFF

function draw_problem(;raff_output=nothing, model_str="logistic")

    datafile = "/tmp/output.txt"

    fp = open(datafile, "r")

    N = parse(Int, readline(fp))

    M = readdlm(fp)

    close(fp)

    x = M[:, 1]
    y = M[:, 2]
    c = M[:, 3]

    true_outliers = findall(c .== 1)

    PyPlot.scatter(x[c .== 0.0], y[c .== 0.0], c=PyPlot.cm."Pastel1"(1.0),
                   marker="o", s=50.0, linewidths=0.2)

    PyPlot.scatter(x[c .== 1.0], y[c .== 1.0], c=PyPlot.cm."Pastel1"(1.0),
                   marker="^", s=50.0, linewidths=0.2, label="Outliers")

    if raff_output != nothing
        
        n, model, modelstr = RAFF.model_list[model_str]

        modl1 = (x) -> model(x, raff_output.solution)

        t = minimum(x):0.01:maximum(x)
        PyPlot.plot(t, modl1.(t), c=PyPlot.cm."Set1"(2.0/9.0))

        # Draw outliers found by RAFF
        
        true_positives = intersect(true_outliers, raff_output.outliers)
        false_positives = setdiff(raff_output.outliers, true_positives)
        
        PyPlot.scatter(x[false_positives], y[false_positives],
                   c=PyPlot.cm."Pastel1"(0.0/9.0), marker="o", s=50.0,
                   linewidths=0.2, label="False positives")
        
        PyPlot.scatter(x[true_positives], y[true_positives],
                       c=PyPlot.cm."Pastel1"(0.0/9.0), marker="^",
                       s=50.0, linewidths=0.2, label="Identified outliers")

    end

    # θSol = [1000.0, 5000.0, 0.2, 3.0]
    # PyPlot.plot(t, modl2.(t), "b--", label="Verdadeiro")

    PyPlot.legend(loc=4)

    PyPlot.show()
    
    PyPlot.savefig("/tmp/figure.png", DPI=100)

end
