using RAFF
using Printf
using Random

using Logging
using Base.CoreLogging

"""

    calc_ratio( model_str::String, np::Int, p::Int,
        sol::Vector{Float64}, ntests::Int=10,
        initguess::Vector{Float64}=nothing, maxms::Int=1,
        fp::IOStream=stdout)

Run a sequence of `ntests` tests for model `model_str` using a dataset
of `np` points and `p` trusted points (i. e. `np - p` outliers). The
points are generated by perturbations around *solution* `sol`. See
function [`generate_noisy_data`](@ref) for more details about the
generation of random points and outliers.

If provided, `initguess` is used as a starting points for
[`raff`](@ref), `maxms` is the number of random initial points used
and `fp` is the output.

"""
function calc_ratio( model_str::String, np::Int, p::Int,
    sol::Vector{Float64}, ntests::Int=10,
    initguess=nothing, maxms::Int=1,fp=stdout; cluster=nothing)
    
    # Set Logging
    global_logger(ConsoleLogger(stdout, Logging.Error))

    n, model, = RAFF.model_list[model_str]

    tmpsol = Vector{Float64}(undef, n)

    large_number = 179424673
    
    # Tests
    tot_tim = 0.0
    n_match = zeros(Int, np + 1)
    n_exact = 0
    n_out = 0
    # Number of correct outliers found
    n_cout = 0
    # Number of false positives found
    n_fp = 0

    for i = 1:ntests

        # Define seed for this run
        Random.seed!(large_number + i)

        if initguess == nothing
            
            x = zeros(Float64, n)
        else

            x = copy(initguess)
            
        end

        tmpsol .= sol

        data, tmpsol, tout = if cluster == nothing
        
            generate_noisy_data(model, n, np, p, tmpsol, (1.0, 30.0))

        else

            generate_clustered_noisy_data(model, n, np, p, tmpsol,
                (1.0, 30.0), cluster)

        end
            
        rsol, t, = @timed praff(model, data[:, 1:end - 1], n;
                                MAXMS=maxms, initguess=x, ftrusted=0.7)

        cnt = 0

        for k in rsol.outliers

            (k in tout) && (cnt += 1)

        end

        n_out += length(rsol.outliers)

        n_cout += cnt

        n_fp += length(rsol.outliers) - cnt

        n_match[cnt + 1] += 1

        (cnt == np - p) && (length(rsol.outliers) == np - p) && (n_exact += 1)

        tot_tim += t

    end

    @printf(fp, "%10s %5d %5d %10.3f %10.3f %10.3f %10.3f %10.2f %5d %5d %5d %10.3f\n",
            model_str, np, p, n_match[np - p + 1] /
            (1.0 * ntests), n_exact / (1.0 * ntests), n_cout / ntests,
            n_fp / ntests, n_out / ntests, n_match[1], n_match[2],
            n_match[3], tot_tim)

    return n_match
    
end


"""

    run_calc_ratio(filename="/tmp/table.txt")

Perform a sequence of tests for different models and a combination of
different number of outliers.

Saves the results in `filename`.

"""
function run_calc_ratio(filename="/tmp/table.txt")

    for (model_str, sol) in [ ("linear", [-200.0, 1000.0]), ("cubic", [0.5, -20.0, 300.0, 1000.0]),
                              ("expon", [5000.0, 4000.0, 0.2]),
                              ("logistic", [6000.0, -5000, -0.2, -3.7]) ]

        for (np, p) in [(10, 9), (10, 8), (100, 99), (100, 90)]

            for maxms in [1, 10, 100, 1000]

                open(filename, "a") do fp

                    calc_ratio(model_str, np, p, sol, 1000, nothing, maxms, fp);

                end

            end

        end

    end

end

"""

    run_calc_ratio_clustered(filename="/tmp/table.txt")

Perform a sequence of tests for different models for solving problems
with clustered outliers.

Saves the results in `filename`.

"""
function run_calc_ratio_clustered(filename="/tmp/table.txt")

    np = 100

    p  = 90

    maxms = 100

    for (model_str, sol) in [ ("linear", [-200.0, 1000.0]), ("cubic", [0.5, -20.0, 300.0, 1000.0]),
                              ("expon", [5000.0, 4000.0, 0.2]),
                              ("logistic", [6000.0, -5000, -0.2, -3.7]) ]

        open(filename, "a") do fp

            calc_ratio(model_str, np, p, sol, 1000, nothing, maxms, fp, cluster=(5.0, 10.0));

        end

    end

end
