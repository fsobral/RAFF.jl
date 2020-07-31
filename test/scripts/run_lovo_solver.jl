using RAFF
using DelimitedFiles
using Printf

# Set Debug for Logging
using Logging
using Base.CoreLogging

function run_lovo_solver(p::Int, initguess=nothing; model_str="logistic", solver=lmlovo)
    
    global_logger(ConsoleLogger(stdout, Logging.Error))

    n, model, modelstr = RAFF.model_list[model_str]
    
    open("/tmp/output.txt") do fp
        
        global N = parse(Int, readline(fp))
        
        global data = readdlm(fp)
        
    end

    if initguess == nothing

        initguess = zeros(Float64, n)

    end

    rsol = solver(model, initguess, data[:, 1:end - 1], n, p)
    
    @printf("Solution found:
            fbest = %f
            p     = %d\n", rsol.f, rsol.p)
    println(rsol.solution)

    return rsol

end
