using RAFF
using DelimitedFiles
using Printf
using Random

using Logging
using Base.CoreLogging

function run_raff(modelStr, np, p, sol, ntests=10, initguess=nothing, maxms=1)
    
    # Set Logging
    global_logger(ConsoleLogger(stdout, Logging.Error))

    n, model, = RAFF.model_list[modelStr]

    tmpsol = Vector{Float64}(undef, n)

    Random.seed!(123456789)

    # Tests
    tot_tim = 0.0
    n_exact = zeros(Int, ntests)
    
    for i = 1:ntests
    
        if initguess == nothing
            
            x = zeros(Float64, n)
        else

            x = copy(initguess)
            
        end

        randn!(tmpsol)

        tmpsol .+= sol
        
        data, tmpsol, tout = generate_noisy_data(model, n, np, p,
                                                 tmpsol, 1.0, 30.0)

        rsol, t, = @timed raff(model, data, n; MAXMS=maxms, initguess=x)

        cnt = 0

        for k in rsol.outliers

            (k in tout) && (cnt += 1)

        end

        n_exact[i] = cnt

        tot_tim += t

    end

    @printf("%10s %5d %5d %10.8f %5d %5d %5d %8.4f\n", modelStr, np, p,
            count(n_exact .== np - p) / (1.0 * ntests),
            count(n_exact .== 0), count(n_exact .== 1),
            count(n_exact .== 2), tot_tim)

    return n_exact
    
end
