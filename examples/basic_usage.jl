using RAFF

#
# This example is the "Basic Usage" example in the documentation
#

# Matrix with data
A=[-2.0  5.0;
   -1.5  3.25;
   -1.0  2.0 ;
   -0.5  1.25;
    0.0  1.0 ;
    0.5  2.55;
    1.0  2.0 ;
    1.5  3.25;
    2.0  5.0 ;]

# Define the model to fit data
model(x, t) = x[1] * t[1]^2 + x[2]

# Number of parameters in the model
n = 2

# Run RAFF
output = raff(model, A, n)

# Print output
println(output)

