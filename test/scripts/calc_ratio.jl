using RAFF
using DelimitedFiles
using Printf
using Random

using Logging
using Base.CoreLogging

function calc_ratio(modelStr, np, p, sol, ntests=10, initguess=nothing, maxms=1, fp=stdout)
    
    # Set Logging
    global_logger(ConsoleLogger(stdout, Logging.Error))

    n, model, = RAFF.model_list[modelStr]

    tmpsol = Vector{Float64}(undef, n)

    Random.seed!(123456789)

    # Tests
    tot_tim = 0.0
    n_match = zeros(Int, ntests)
    n_exact = 0
    n_out = 0
    
    for i = 1:ntests
    
        if initguess == nothing
            
            x = zeros(Float64, n)
        else

            x = copy(initguess)
            
        end

        randn!(tmpsol)

        tmpsol .+= sol
        
        data, tmpsol, tout = generate_noisy_data(model, n, np, p,
                                                 tmpsol, (1.0, 30.0))

        rsol, t, = @timed raff(model, data, n; MAXMS=maxms, initguess=x)

        cnt = 0

        for k in rsol.outliers

            (k in tout) && (cnt += 1)

        end

        n_out += length(rsol.outliers)

        n_match[i] = cnt

        (cnt == np - p) && (length(rsol.outliers) == np - p) && (n_exact += 1)

        tot_tim += t

    end

    @printf(fp, "%10s %5d %5d %10.8f %10.8f %10.2f %5d %5d %5d %8.4f\n", modelStr, np, p,
            count(n_match .== np - p) / (1.0 * ntests),
            n_exact / (1.0 * ntests),
            n_out / ntests,
            count(n_match .== 0), count(n_match .== 1),
            count(n_match .== 2), tot_tim)

    return n_match
    
end

function run_calc_ratio()

    for (modelStr, sol) in [("linear", [1.0, -7.0])] #, "cubic", "expon", "logistic"]

        for (np, p) in [(100, 99), (100, 90)]
            #(10, 9), (10, 8) , (100, 99), (100, 90), (1000, 1)]

            for maxms in [1, 10, 100, 1000]

                open("/tmp/table.txt", "a") do fp

                    calc_ratio(modelStr, np, p, sol, 100, nothing, maxms, fp);

                end

            end

        end

    end

end
