using DelimitedFiles
using PyPlot

datafile = "/tmp/output.txt"

fp = open(datafile, "r")

N = parse(Int, readline(fp))

M = readdlm(fp)

close(fp)

x = M[:, 1]
y = M[:, 2]
c = M[:, 3]

sol = [699.522, 5476.89, 0.16228, 14.3365]

#model = (x, t) -> x[1] .* t.^3 + x[2] .* t.^2 + x[3] .* t .+ x[4]
#model = (x, t) -> x[1] .+ x[2] .* exp.(- x[3] .* t)
model = (x, t) -> x[1] .+ x[2] ./ (1.0 .+ exp.(- x[3] .* t .- x[4]))

t = minimum(x):0.01:maximum(x)
PyPlot.plot(t, model(sol, t), "r")

xSol = [1000.0, 5000.0, 0.2, 15.0]
PyPlot.plot(t, model(xSol, t), "b--")


PyPlot.scatter(x, y, c=c, marker="o", s=50.0, linewidths=0.2,
               cmap=PyPlot.cm["Paired"], alpha=0.6)

PyPlot.show()
