using RAFF
using DelimitedFiles
using PyPlot

#
# This example shows how solve a problem given by data in a
# file. Also, we show how to use the `model_list` utility structure to
# retrieve pre-defined models in RAFF.
#

# Retrieve the number of parameters, the model function and the model
# function string (not used) for the "Cubic" model.
n, model, modelstr = RAFF.model_list["cubic"]

# Read data from file
open("C1.txt") do fp
        
    global data = readdlm(fp)
        
end

# Set starting point
initguess = [0.0, 0.0, 0.0, 0.0]

# Set the number of multistart iterations.
maxms = 10

# Call RAFF. In this case, the model is not a multivariate one. The
# last column of the file indicates which observation is the outlier.
rsol = raff(model, data[:, 1:end - 1], n; MAXMS=maxms, initguess=initguess)
    
println(rsol)

# Now we plot the solution of the problem

x = data[:, 1]
y = data[:, 2]
c = data[:, 3]

# Plot the obtained solution
modl1 = (t) -> model(rsol.solution, t)

t = minimum(x):0.01:maximum(x)
PyPlot.plot(t, modl1.(t), "r", label="RAFF")

# Plot the 'true' solutios, i. e., the one use to generate the
# problem.
xSol = [2.0, 0, -4.0, -10]
modl2 = (t) -> model(xSol, t)

PyPlot.plot(t, modl2.(t), "b--", label="True Solution")

PyPlot.legend(loc=4)

# Color the obtained outliers
c[rsol.outliers] .= 5.0

# Plot all the points
PyPlot.scatter(x, y, c=c, marker="o", s=50.0, linewidths=0.2,
               cmap=PyPlot.cm."Paired", alpha=0.6)

PyPlot.show()
