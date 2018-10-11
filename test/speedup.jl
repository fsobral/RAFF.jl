using RAFF
using Random
using Distributed
using Printf
using ArgParse

# Set Debug for Logging
using Logging
using Base.CoreLogging

global_logger(ConsoleLogger(stdout, Logging.Error))

s = ArgParseSettings()

@add_arg_table s begin

    "type"
    help = "Sequential(s) or parallel(p)"
    arg_type = String
    required = true

    "np"
    arg_type = Int
    help = "Number of points"
    required = true

    "nw"
    arg_type = Int
    help = "Number of workers"

end

parsed_args = parse_args(ARGS, s)

n = 2

np = parsed_args["np"]
    
p = Int(0.7 * np)
    
answer = [2.0, -0.5]

if parsed_args["type"] == "p"

    ###################
    # Distributed run #
    ###################

    nw = parsed_args["nw"]

    addprocs(nw)
    
    @everywhere using RAFF

    # Set Debug for Logging
    @everywhere using Logging
    @everywhere using Base.CoreLogging

    @everywhere global_logger(ConsoleLogger(stdout, Logging.Error))

    @everywhere gmodel!(x, t_, g) = begin
        g[1] = exp(t_ * x[2])
        g[2] = t_ * x[1] * exp(t_ * x[2]);
    end

    @everywhere model(x, t) = x[1] * exp(t * x[2])

    # First run
    praff(model, gmodel!, zeros(1, 3), n)
    
    # Set the seed for generating the same data
    Random.seed!(123456789)

    data, xSol = RAFF.generateNoisyData(model, n, np, p, answer)

    val, t, bytes, gctime, memallocs = try

        @timed praff(model, gmodel!, data, n)

    catch e

        ([1.0e+99, 1.0e+99], 1.0e+99, -1), -1.0, 0, -1.0, nothing

    end

    sol, fsol, psol = val
        
    @printf("%2d %10d %10.5s %15.8e %15.8e %15.8e %6d\n", nw, np,
            t, sol[1], sol[2], fsol, psol)

    rmprocs(workers())

else

    ##############
    # Serial run #
    ##############

    gmodel!(x, t_, g) = begin
        g[1] = exp(t_ * x[2])
        g[2] = t_ * x[1] * exp(t_ * x[2]);
    end

    model(x, t) = x[1] * exp(t * x[2])

    # First run
    raff(model, gmodel!, zeros(1, 3), n)

    # Set the seed for generating the same data
    Random.seed!(123456789)
    
    data, xSol = RAFF.generateNoisyData(model, n, np, p, answer)
    
    val, t, bytes, gctime, memallocs = try

        @timed raff(model, gmodel!, data, n)

    catch e

        ([1.0e+99, 1.0e+99], 1.0e+99, -1), -1.0, 0, -1.0, nothing

    end

    sol, fsol, psol = val
    
    @printf("%2d %10d %10.5s %15.8e %15.8e %15.8e %6d\n", 0, np,
            t, sol[1], sol[2], fsol, psol)

end
