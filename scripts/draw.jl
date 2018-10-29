using DelimitedFiles
using PyPlot
using ArgParse

datafile = "/tmp/L5.txt"

M = readdlm(datafile)

x = M[:, 1]
y = M[:, 2]
c = [Int(i) for i in (M[:, 3] .!= 0.0)] * 1.0

PyPlot.scatter(x, y, c=c, marker="o", s=50.0, linewidths=0.2,
               cmap=PyPlot.cm["Paired"], alpha=0.6)

PyPlot.show()
