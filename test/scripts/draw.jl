using DelimitedFiles
using PyPlot
using RAFF

datafile = "/tmp/output.txt"

fp = open(datafile, "r")

N = parse(Int, readline(fp))

M = readdlm(fp)

close(fp)

x = M[:, 1]
y = M[:, 2]
c = M[:, 3]

# sol = [699.522, 5476.89, 0.16228, 2.32653]
sol = [533.077, 5389.48, 0.153425, 2.04606]

n, model, modelstr = RAFF.model_list["logistic"]

modl1 = (x) -> model(x, sol)
modl2 = (x) -> model(x, θSol)

# t = minimum(x):0.01:maximum(x)
# PyPlot.plot(t, modl1.(t), "r", label="Ajustado RAFF")

# θSol = [1000.0, 5000.0, 0.2, 3.0]
# PyPlot.plot(t, modl2.(t), "b--", label="Verdadeiro")

# PyPlot.legend(loc=4)

PyPlot.scatter(x, y, c=c, marker="o", s=50.0, linewidths=0.2,
               cmap=PyPlot.cm."Paired", alpha=0.6)



PyPlot.show()

PyPlot.savefig("/tmp/figure.png", DPI=100)
