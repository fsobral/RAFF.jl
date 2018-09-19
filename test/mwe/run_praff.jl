# This is a minimal working example for running praff on several
# workers in Julia and saving them into a shared matrix

@everywhere using RAFF
@everywhere using Base.CoreLogging

@everywhere CoreLogging.disable_logging(CoreLogging.Info)

@everywhere gmodel!(x, t_, g) = begin
    g[1] = exp(t_ * x[2])
    g[2] = t_ * x[1] * exp(t_ * x[2]);
end

@everywhere model(x, t) = x[1] * exp(t * x[2])

n = 2

np = 10

p = 7

data, xSol = RAFF.generateNoisyData(model, n, np, p)

praff(model, gmodel!, data, n)
