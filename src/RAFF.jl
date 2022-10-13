"""

`RAFF.jl` is a Jula package.

"""
module RAFF

# Dependencies
using Distributed
using ForwardDiff
using LinearAlgebra
using Statistics
using Printf
using Random
using SharedArrays
using Logging

export raff, praff

const raff_logger = Ref{AbstractLogger}()
const lm_logger = Ref{AbstractLogger}()

function __init__()

    # Set RAFF logger
    raff_logger.x = ConsoleLogger(stdout, Logging.Error)
    lm_logger.x = ConsoleLogger(stdout, Logging.Error)
    
end

# Load code

include("raffoutput.jl")

include("utils.jl")

include("dutils.jl")

include("generator.jl")

include("inner_solvers.jl")

include("smartqr.jl")

"""

    raff(model::Function, data::Array{Float64, 2}, n::Int; kwargs...)

    raff(model::Function, gmodel!::Function, data::Array{Float64, 2},
        n::Int; MAXMS::Int=1, SEEDMS::Int=123456789,
        initguess::Vector{Float64}=zeros(Float64, n),
        noutliers::Int=-1, ftrusted::Union{Float64,
        Tuple{Float64, Float64}}=0.5,
        inner_solver::Function=lmlovo, inner_solver_params...)

Robust Algebric Fitting Function (RAFF) algorithm. This function uses
a voting system to automatically find the number of trusted data
points to fit the `model`.

  - `model`: function to fit data. Its signature should be given by

        model(x, θ)

    where `x` is the multidimensional argument and `θ` is the
    `n`-dimensional vector of parameters

  - `gmodel!`: gradient of the model function. Its signature should be
    given by

        gmodel!(g, x, θ)

    where `x` is the multidimensional argument, `θ` is the
    `n`-dimensional vector of parameters and the gradient is written
    in `g`.

  - `data`: data to be fit. This matrix should be in the form

        x11 x12 ... x1N y1
        x21 x22 ... x2N y2
        :

    where `N` is the dimension of the argument of the model
    (i.e. dimension of `x`).

  - `n`: dimension of the parameter vector in the model function

The optional arguments are

  - `MAXMS`: number of multistart points to be used
  - `SEEDMS`: integer seed for random multistart points
  - `initialguess`: a good guess for the starting point and for
    generating random points in the multistart strategy
  - `noutliers`: integer describing the maximum expected number of
    outliers. The default is *half*. *Deprecated*.
  - `ftrusted`: float describing the minimum expected percentage of
    trusted points. The default is *half* (0.5). Can also be a
    Tuple of the form `(fmin, fmax)` percentages of trusted points.
  - `inner_solver`: solver to be used for the least square problems.
    By default, uses [`lmlovo`](@ref). This function has the following
    mandatory parameters

        inner_solver(model, gmodel!, θ, data, n, p;
                     inner_solver_params...) = RAFFOutput

  - `inner_solver_params...`: the remaining parameters will be sent
    as optional arguments to the `inner_solver`

Returns a [`RAFFOutput`](@ref) object with the best parameter found.

"""
function raff(model::Function, gmodel!::Function,
    data::Array{Float64, 2}, n::Int; MAXMS::Int=1,
    SEEDMS::Int=123456789, initguess::Vector{Float64}=zeros(Float64,
    n), noutliers::Int=-1,
    ftrusted::Union{Float64, Tuple{Float64, Float64}}=0.5,
    inner_solver::Function=lmlovo, inner_solver_params...)

    np, = size(data)
    
    if noutliers != -1

        Base.depwarn("Optional argument `noutliers::Int` is deprecated and will be removed from future releases. Use `ftrusted::Float` or `itrusted::Tuple{Float, Float}` instead.", :raff)

        ftrusted = (noutliers >= 0) ? max(0, np - noutliers) / np : 0.5

    end
    
    # Initialize random generator
    seedMS = MersenneTwister(SEEDMS)

    # Define interval for trusted points

    pliminf, plimsup = try

        check_ftrusted(ftrusted, np)

    catch e

        with_logger(raff_logger.x) do

            @error("Error in optional parameter `ftrusted`.", e)

        end

        return RAFFOutput()

    end
    
    lv = plimsup - pliminf + 1
    
    sols = Vector{RAFFOutput}(undef, lv)

    for i = pliminf:plimsup

        vbest = RAFFOutput(initguess, i)
        
        ind = i - pliminf + 1
        
        for j = 1:MAXMS

            with_logger(raff_logger.x) do
                
                @debug("Running $(string(inner_solver)) for p = $(i). Repetition $(j).")
                
            end
            
            # Starting point
            θ = randn(seedMS, Float64, n)
            θ .= θ .+ initguess
        
            # Call function and store results
            sols[ind] = inner_solver(model, gmodel!, θ, data, n, i;
                                     inner_solver_params...)

            # Update the best point and functional value
            (sols[ind].status == 1) && (sols[ind].f < vbest.f) &&
                (vbest = sols[ind])
            
        end

        sols[ind] = vbest

        with_logger(raff_logger.x) do
                
            @debug("Best solution for p = $(i).", vbest.solution)
                
        end
        

    end

    # Count the total number of iterations, and function and Jacobian
    # evaluations.

    nf = 0
    nj = 0
    ni = 0

    for s in sols

        nf += s.nf
        nj += s.nj
        ni += s.iter

    end

    # Apply the filter and the voting strategy to all the solutions
    # found

    votsis = with_logger(raff_logger.x) do

        voting_strategy(model, data, sols, pliminf, plimsup)

    end

    mainind = findlast(x->x == maximum(votsis), votsis)

    s = sols[mainind]
    
    return RAFFOutput(s.status, s.solution, ni, s.p, s.f, nf, nj, s.outliers)
    
end

function raff(model::Function, data::Array{Float64, 2}, n::Int; kwargs...)

    # Define closures for derivative and initializations

    # 'x' is considered as global for this function
    model_cl(θ) = model(x, θ)
    
    grad_model!(g, x_, θ) = begin

        global x = x_
        
        return ForwardDiff.gradient!(g, model_cl, θ)

    end

    return raff(model, grad_model!, data, n; kwargs...)

end

"""

    praff(model::Function, data::Array{Float64, 2}, n::Int; kwargs...)

    praff(model::Function, gmodel!::Function, data::Array{Float64, 2},
        n::Int; MAXMS::Int=1, SEEDMS::Int=123456789, batches::Int=1,
        initguess::Vector{Float64}=zeros(Float64, n),
        noutliers::Int=-1, ftrusted::Union{Float64,
        Tuple{Float64, Float64}}=0.5,
        inner_solver::Function=lmlovo, inner_solver_params...)

Multicore distributed version of RAFF. See the description of the
[`raff`](@ref) function for the main (non-optional) arguments. All the
communication is performed by channels.

This function uses all available **local** workers to run RAFF
algorithm. Note that this function does not use *Tasks*, so all the
parallelism is based on the [Distributed](https://docs.julialang.org/en/latest/manual/parallel-computing/#Multi-Core-or-Distributed-Processing-1) package.

The optional arguments are

  - `MAXMS`: number of multistart points to be used
  - `SEEDMS`: integer seed for random multistart points
  - `batches`: size of batches to be send to each worker
  - `initguess`: starting point to be used in the multistart procedure
  - `noutliers`: integer describing the maximum expected number of
    outliers. The default is *half*. *Deprecated*.
  - `ftrusted`: float describing the minimum expected percentage of
    trusted points. The default is *half* (0.5). Can also be a
    Tuple of the form `(fmin, fmax)` percentages of trusted points.
  - `inner_solver`: solver to be used for the least square problems.
    By default, uses [`lmlovo`](@ref). This function has the following
    mandatory parameters

        inner_solver(model, gmodel!, θ, data, n, p;
                     inner_solver_params...) = RAFFOutput

  - `inner_solver_params...`: the remaining parameters will be sent
    as optional arguments to the `inner_solver`

Returns a [`RAFFOutput`](@ref) object containing the solution.

"""
function praff(model::Function, gmodel!::Function,
               data::Array{Float64, 2}, n::Int; MAXMS::Int=1,
               SEEDMS::Int=123456789, batches::Int=1,
               initguess::Vector{Float64}=zeros(Float64, n),
               noutliers::Int=-1,
               ftrusted::Union{Float64, Tuple{Float64, Float64}}=0.5,
               inner_solver::Function=lmlovo, inner_solver_params...)

    np, = size(data)
    
    if noutliers != -1

        Base.depwarn("Optional argument `noutliers::Int` is deprecated and will be removed from future releases. Use `ftrusted::Float` or `itrusted::Tuple{Float, Float}` instead.", :raff)

        ftrusted = (noutliers >= 0) ? max(0, np - noutliers) / np : 0.5

    end
    
    # Initialize random generator
    seedMS = MersenneTwister(SEEDMS)

    # Define interval for trusted points

    pliminf, plimsup = try

        check_ftrusted(ftrusted, np)

    catch e

        with_logger(raff_logger.x) do

            @error("Error in optional parameter `ftrusted`.", e)

        end

        return RAFFOutput()

    end

    lv = plimsup - pliminf + 1
    
    # Create a RemoteChannel to receive solutions
    bqueue = RemoteChannel(() -> Channel{Vector{Float64}}(div(lv, 2)))
    # TODO: Check a smart way for not creating a large channel
    squeue = RemoteChannel(() -> Channel{RAFFOutput}(lv))
    # Create another channel to assign tasks
    tqueue = RemoteChannel(() -> Channel{UnitRange{Int}}(0))

    # This command selects only nodes which are local to myid()
    curr_workers = workers()

    futures = Vector{Future}(undef, length(curr_workers))

    # Start updater Task
    # This task is not needed up to now.
    # @async with_logger(()-> update_best(bqueue, bestx), raff_logger)
        
    # Start workers Tasks (CPU intensive)
    with_logger(raff_logger.x) do

        @debug("Workers", curr_workers)

    end
    
    for (i, t) in enumerate(curr_workers)

        futures[i] = @spawnat(t, with_logger( ()-> try

            @debug("Creating worker $(t).")
                                              
            consume_tqueue(bqueue, tqueue, squeue,
                           model, gmodel!, data, n, pliminf,
                           plimsup, MAXMS, seedMS, initguess,
                           inner_solver; inner_solver_params...)
            catch e
                                              
               @error("Unable to start worker $(t).", e)

            end, raff_logger.x
        ))

        with_logger(raff_logger.x) do

            @debug("Created worker $(t).")

        end
    end

    # Check asynchronously if there is at least one live worker
    @async with_logger(
        () -> check_and_close(bqueue, tqueue, squeue, futures),
        raff_logger.x)

    # Populate the task queue with jobs
    
    for p = pliminf:batches:plimsup

        try
            
            put!(tqueue, p:min(plimsup, p + batches - 1))

        catch e

            with_logger(raff_logger.x) do

                @warn("Tasks queue prematurely closed while inserting tasks. Will exit.")

            end

            break

        end

        with_logger(raff_logger.x) do

            @debug("Added problem $(p) to tasks queue.")

        end
                

    end

    # The task queue can be closed, since all the problems have been
    # read, due to the size 0 of this channel
    close(tqueue)

    with_logger(raff_logger.x) do
        
        @debug("Waiting for workers to finish.")

    end

    # Create a vector of solutions to store the results from workers
    sols = Vector{RAFFOutput}(undef, lv)

    for i in 1:lv

        try

            rout = take!(squeue)

            sols[rout.p - pliminf + 1] = rout

            with_logger(raff_logger.x) do

                @debug("Stored solution for p=$(rout.p).")

            end

        catch e

            with_logger(raff_logger.x) do

                @error("Error when retrieving solutions.", e)

            end
            
        end

    end
    
    close(bqueue)

    close(squeue)

    # Count the total number of iterations, and function and Jacobian
    # evaluations.

    nf = 0
    nj = 0
    ni = 0

    for s in sols

        nf += s.nf
        nj += s.nj
        ni += s.iter

    end

    # Apply the filter and the voting strategy to all the solutions
    # found

    votsis = with_logger(raff_logger.x) do

        voting_strategy(model, data, sols, pliminf, plimsup)

    end

    mainind = findlast(x->x == maximum(votsis), votsis)
    
    s = sols[mainind]

    return RAFFOutput(s.status, s.solution, ni, s.p, s.f, nf, nj, s.outliers)

end

function praff(model::Function, data::Array{Float64, 2}, n::Int; kwargs...)

    # Define closures for derivative and initializations

    # 'x' is considered as global parameter for this function
    model_cl(θ) = model(x, θ)
    
    grad_model!(g, x_, θ) = begin

        global x = x_
        
        return ForwardDiff.gradient!(g, model_cl, θ)

    end

    return praff(model, grad_model!, data, n; kwargs...)

end

end
