using RAFF
using Distributed

#
# This example is the "Parallel running" example from the
# documentation.
#

# Add 3 worker processes
addprocs(3)

# Distribute RAFF, the model and its gradien (the gradient is not
# mandatory, just an example)
@everywhere using RAFF

@everywhere function model(x, t)
   x[1] * t[1]^2 + x[2]
end

@everywhere function gmodel!(x, t, g)
    g[1] = t[1]^2
    g[2] = 1.0
end

# Define the data matrix
A=[-2.0  5.0;
  -1.5  3.25;
  -1.0  2.0 ;
  -0.5  1.25;
   0.0  1.0 ;
   0.5  2.55;
   1.0  2.0 ;
   1.5  3.25;
   2.0  5.0 ;];

# Number of parameters of the model
n = 2

# Run Parallel RAFF
output = praff(model, gmodel!, A, n)

# Print output
println(output)
