# This is a minimal working example for running praff on several
# workers in Julia and saving them into a shared matrix

using Random

@everywhere using RAFF
@everywhere using Base.CoreLogging

@everywhere CoreLogging._min_enabled_level[] = CoreLogging.Warn + 1

@everywhere gmodel!(x, t_, g) = begin
    g[1] = exp(t_ * x[2])
    g[2] = t_ * x[1] * exp(t_ * x[2]);
end

@everywhere model(x, t) = x[1] * exp(t * x[2])

n = 2

np = 100

p = 70

# Set the seed for generating the same data
Random.seed!(123456789)

data, xSol = RAFF.generateNoisyData(model, n, np, p)

praff(model, gmodel!, data, n, MAXMS=2)

println("True sol: $(xSol).")
