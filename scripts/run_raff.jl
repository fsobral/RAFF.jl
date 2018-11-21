using RAFF
using DelimitedFiles
using Printf

# Set Debug for Logging
using Logging
using Base.CoreLogging

function run_raff(maxms=1, initguess=nothing)
    
    global_logger(ConsoleLogger(stdout, Logging.Debug))

    # model = (x, t) -> x[1] * t^3 + x[2] * t^2 + x[3] * t + x[4]
    # model = (x, t) -> x[1] + x[2] * exp(- x[3] * t)
    model = (x, t) -> x[1] + x[2] / (1.0 + exp(- x[3] * (t - x[4])))

    n = 4
    
    open("/tmp/output.txt") do fp
        
        global N = parse(Int, readline(fp))
        
        global data = readdlm(fp)
        
    end

    if initguess == nothing

        initguess = zeros(Float64, n)

    end

    xbest, fbest, p = raff(model, data, n; MAXMS=maxms,
                           initguess=initguess)
    
    @printf("Solution found:
            fbest = %f
            p     = %d\n", fbest, p)
    println(xbest)

end
