using DelimitedFiles
using PyPlot
using ArgParse

datafile = "/tmp/output.txt"

M = readdlm(datafile)

x = M[:, 1]
y = M[:, 2]
c = [Int(i) for i in (M[:, 3] .!= 0.0)] * 1.0

sol = [2.03607, 0.7784, -7.0366, -68.5755]
model = (x, t) -> x[1] .* t.^3 + x[2] .* t.^2 + x[3] .* t .+ x[4]
t = minimum(x):0.01:maximum(x)
PyPlot.plot(t, model(sol, t), "r")

#xSol = [2.0, 0.0, -4.0, -10.0]
#PyPlot.plot(t, model(xSol, t), "b")

PyPlot.scatter(x, y, c=c, marker="o", s=50.0, linewidths=0.2,
               cmap=PyPlot.cm["Paired"], alpha=0.6)

PyPlot.show()
