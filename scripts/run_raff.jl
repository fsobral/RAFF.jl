using RAFF
using DelimitedFiles
using Printf

# Set Debug for Logging
using Logging
using Base.CoreLogging

function run_raff(maxms=1, initguess=nothing)
    
    global_logger(ConsoleLogger(stdout, Logging.Error))

    # model = (x, t) -> x[1] * t^3 + x[2] * t^2 + x[3] * t + x[4]
    # model = (x, t) -> x[1] + x[2] * exp(- x[3] * t)
    model = (x, t) -> x[1] + x[2] / (1.0 + exp(- x[3] * t[1] - x[4]))

    n = 4
    
    open("/tmp/output.txt") do fp
        
        global N = parse(Int, readline(fp))
        
        global data = readdlm(fp)
        
    end

    if initguess == nothing

        initguess = zeros(Float64, n)

    end

    rsol = raff(model, data[:, 1:2], n; MAXMS=maxms, initguess=initguess, ε=1.0e-6)
    
    @printf("Solution found:
            fbest = %f
            p     = %d\n", rsol.f, rsol.p)
    println(rsol.solution)

end
