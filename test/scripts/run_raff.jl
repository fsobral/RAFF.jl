using RAFF
using DelimitedFiles
using Printf

# Set Debug for Logging
using Logging
using Base.CoreLogging

function run_raff(maxms=1, initguess=nothing;
                  model_str="logistic", foutliers=0.5)
    
    n, model, modelstr = RAFF.model_list[model_str]

    open("/tmp/output.txt") do fp
        
        global N = parse(Int, readline(fp))
        
        global data = readdlm(fp)
        
    end

    if initguess == nothing

        initguess = zeros(Float64, n)

    end

    rsol = raff(model, data[:, 1:end - 1], n; MAXMS=maxms, initguess=initguess,
                noutliers=Int(round(foutliers * size(data)[1])))
    
    @printf("Solution found:
            fbest = %f
            p     = %d\n", rsol.f, rsol.p)
    println(rsol.solution)

    return rsol

end
