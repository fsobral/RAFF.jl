using RAFF

#
# This example is the "Multivariate" example in the documentation. It
# also shows some parameters of RAFF. To check all the possible
# parameters, please refer to the documentation.
#

# Matrix with data. Now the experiments has two variables
data = [1.0 1.0    2.0
        0.0 0.0    4.0
        7.0 1.5   -4.5
        2.0 2.0  -17.0 # outlier
        0.0 8.6   -4.6]

# Define the model to fit data
model(x, θ) = θ[1] * x[1] + θ[2] * x[2] + θ[3]            

# Number of parameters in the model (dimension of θ)
n = 3

# Run RAFF. Uncomment the other piece of code, in order to run RAFF
# with different options.

# output = raff(model, data, 3; MAXMS=1, initguess=[0.5,0.5,0.5], ε=1.0e-10)
output = raff(model, data, n)

# Print output
println(output)
