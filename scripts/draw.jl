using DelimitedFiles
using PyPlot
using ArgParse

datafile = "/tmp/output.txt"

fp = open(datafile, "r")

N = parse(Int, readline(fp))

M = readdlm(fp)

close(fp)

x = M[:, 1]
y = M[:, 2]
c = M[:, 3]

#sol = [2.07816, 1.08284, -5.55829, -20.9451]

#model = (x, t) -> x[1] .* t.^3 + x[2] .* t.^2 + x[3] .* t .+ x[4]
model = (x, t) -> x[1] .* exp.(- x[2] .* t)

t = minimum(x):0.01:maximum(x)
#PyPlot.plot(t, model(sol, t), "r")

xSol = [2.0, 2.0]
PyPlot.plot(t, model(xSol, t), "b--")


PyPlot.scatter(x, y, c=c, marker="o", s=50.0, linewidths=0.2,
               cmap=PyPlot.cm["Paired"], alpha=0.6)

PyPlot.show()
