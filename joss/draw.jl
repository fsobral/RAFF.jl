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

sol = [779.616, 5315.36, 0.174958, 2.52718]

n, model, modelstr = RAFF.model_list["logistic"]

modl1 = (x) -> model(x, sol)
modl2 = (x) -> model(x, θSol)

t = minimum(x):0.01:maximum(x)
PyPlot.plot(t, modl1.(t), "r", label="RAFF")

θSol = [1000.0, 5000.0, 0.2, 3.0]
PyPlot.plot(t, modl2.(t), "b--", label="True")

θLS = [-959.07, 8151.03, 0.0927191, 0.940711]
modl3 = (x) -> model(x, θLS)

PyPlot.plot(t, modl3.(t), "g-.", label="Least Squares")

PyPlot.legend(loc=4)

PyPlot.scatter(x[c .== 0.0], y[c .== 0.0], color=PyPlot.cm."Paired"(0.1), marker="o", s=50.0, linewidths=0.2,
               alpha=1.0)

PyPlot.scatter(x[c .== 1.0], y[c .== 1.0], color=PyPlot.cm."Paired"(0.5), marker="^", s=50.0, linewidths=0.2,
               alpha=1.0)

PyPlot.show()

PyPlot.savefig("/tmp/figure.png", DPI=100)
