using RAFF
using Random
using Distributed
using Printf

# Set Debug for Logging
using Logging
using Base.CoreLogging

global_logger(ConsoleLogger(stdout, Logging.Error))

n = 2

answer = [2.0, -0.5]

###################
# Distributed run #
###################

for nw = 0:Sys.CPU_THREADS

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
    
    for i = 1:3

        np = 10^i

        p = Int(0.7 * np)

        # Set the seed for generating the same data
        Random.seed!(123456789)

        data, xSol = RAFF.generateNoisyData(model, n, np, p, answer)

        val, t, bytes, gctime, memallocs = @timed praff(model, gmodel!, data, n)

        sol, fsol, psol = val
        
        @printf("%2d 10^%d %10.5s %15.8e %15.8e %15.8e %6d\n", nw, i,
                t, sol[1], sol[2], fsol, psol)

    end

    rmprocs(workers())

end

@printf("\n")

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

for i = 1:3
    
    np = 10^i
    
    p = Int(0.7 * np)
    
    # Set the seed for generating the same data
    Random.seed!(123456789)
    
    data, xSol = RAFF.generateNoisyData(model, n, np, p, answer)
    
    val, t, bytes, gctime, memallocs = @timed raff(model, gmodel!, data, n)

    sol, fsol, psol = val
        
    @printf("%2d 10^%d %10.5s %15.8e %15.8e %15.8e %6d\n", 0, i,
            t, sol[1], sol[2], fsol, psol)

end
