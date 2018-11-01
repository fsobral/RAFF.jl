using RAFF
using DelimitedFiles
using Printf

# Set Debug for Logging
using Logging
using Base.CoreLogging

global_logger(ConsoleLogger(stdout, Logging.Debug))

model(x, t) = x[1] * t^3 + x[2] * t^2 + x[3] * t + x[4]

n = 4

np = 10

p = 9

data = readdlm("/tmp/output.txt")

xbest, fbest, p = raff(model, data, n)

@printf("Solution found:
        fbest = %f
        p     = %d\n", fbest, p)
println(xbest)
