# This is a minimal working example for running LMlovo on several
# workers in Julia and saving them into a shared matrix

using Revise
using Distributed
using SharedArrays

@everywhere using RAFF
@everywhere using Base.CoreLogging

@everywhere CoreLogging._min_enabled_level[] = CoreLogging.Warn + 1

@everywhere gmodel!(x, t_, g) = begin
    g[1] = exp(t_ * x[2])
    g[2] = t_ * x[1] * exp(t_ * x[2]);
end

@everywhere model(x, t) = x[1] * exp(t * x[2])

n = 2

np = 1000

p = 700

data, xSol = RAFF.generateNoisyData(model, n, np, p)

result = SharedArray{Float64, 2}(n, np)

f = @sync @distributed for i = 1:np

    # Starting point
    x = zeros(Float64, 2)

    # Call function and store results
    s, x, iter, p, f = LMlovo(model, gmodel!, x, data, 2, i, MAXITER=1000)

    result[:, i] .= x

    println("Finished. p = $(p) and f = $(f).")
    
end
