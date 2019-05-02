# This is a minimal working example for running lmlovo function

using Random
using RAFF

# Set Logging.Debug for Logging
using Logging
using Base.CoreLogging

global_logger(ConsoleLogger(stdout, Logging.Debug))

gmodel!(x, t_, g) = begin
    g[1] = exp(t_ * x[2])
    g[2] = t_ * x[1] * exp(t_ * x[2]);
end

model(x, t) = x[1] * exp(t * x[2])

n = 2

np = 100

p = 70

# Set the seed for generating the same data
Random.seed!(123456789)

data, xSol = RAFF.generateNoisyData(model, n, np, p, [2.0, -0.5])

@time status, x, iter, p_, f = lmlovo(model, data, n, p)

println("True sol: $(xSol).")
println("Found: $x.")
