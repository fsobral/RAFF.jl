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

export lmlovo, raff, praff

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

"""
    lmlovo(model::Function [, θ::Vector{Float64} = zeros(n)], data::Array{Float64, 2},
           n::Int, p::Int [; kwargs...])

    lmlovo(model::Function, gmodel!::Function [, θ::Vector{Float64} = zeros(n)],
           data::Array{Float64,2}, n::Int, p::Int [; MAXITER::Int=200,
           ε::Float64=10.0^-4])

Fit the `n`-parameter model `model` to the data given by matrix
`data`. The strategy is based on the LOVO function, which means that
only `p` (0 < `p` <= rows of `data`) points are trusted. The
Levenberg-Marquardt algorithm is implemented in this version.

Matriz `data` is the data to be fit. This matrix should be in the form

    x11 x12 ... x1N y1
    x21 x22 ... x2N y2
    :

where `N` is the dimension of the argument of the model
(i.e. dimension of `x`).

If `θ` is provided, then it is used as the starting point.

The signature of function `model` should be given by

    model(x::Union{Vector{Float64}, SubArray}, θ::Vector{Float64})

where `x` are the variables and `θ` is a `n`-dimensional vector of
parameters. If the gradient of the model `gmodel!`

    gmodel! = (g::SubArray, x::Union{Vector{Float64}, SubArray},
               θ::Vector{Float64})

is not provided, then the function ForwardDiff.gradient! is called to
compute it.  **Note** that this choice has an impact in the
computational performance of the algorithm. In addition, if
`ForwardDiff.jl` is being used, then one **MUST** remove the signature
of vector `θ` from function `model`.

The optional arguments are

  - `MAXITER`: maximum number of iterations
  - `ε`: tolerance for the gradient of the function

Returns a [`RAFFOutput`](@ref) object.

"""
function lmlovo(model::Function, gmodel!::Function, θ::Vector{Float64},
                data::Array{Float64,2}, n::Int, p::Int;
                MAXITER::Int=200, ε::Float64=10.0^-4)

    @assert(n > 0, "Dimension should be positive.")
    @assert(p >= 0, "Trusted points should be nonnegative.")
    
    npun, = size(data)

    with_logger(lm_logger.x) do
        
        @debug("Size of data matrix ", size(data))

    end

    # Counters for calls to F and its Jacobian
    nf = 0

    nj = 0
    
    (p == 0) && return RAFFOutput(1, θ, 0, p, 0.0, nf, nj, [1:npun;])

    # Main function - the LOVO function
    LovoFun = let

        npun_::Int = npun
        
        ind::Vector{Int} = Vector{Int}(undef, npun_)
        
        F::Vector{Float64} = Vector{Float64}(undef, npun_)
        
        p_::Int = p
        
        # Return a ordered set index and lovo value
        (θ) -> begin

            nf += 1
            
            @views for i = 1:npun_
                F[i] = (model(data[i,1:(end - 1)], θ) - data[i, end])^2
            end
            
            indF, orderedF = sort_fun!(F, ind, p_)
            
            return indF, sum(orderedF)
        end
        
    end
    
    # Residue and Jacobian of residue
    val_res::Vector{Float64} = Vector{Float64}(undef, p)
    
    jac_res::Array{Float64, 2} = Array{Float64}(undef, p, n)
    
    # This function returns the residue and Jacobian of residue
    ResFun!(θ::Vector{Float64}, ind, r::Vector{Float64},
            rJ::Array{Float64, 2}) = begin

       nj += 1
                
       for (k, i) in enumerate(ind)
            
            x = data[i, 1:(end - 1)]
            
            r[k] = model(x, θ) - data[i, end]
            
            v = @view(rJ[k, :])
            
            gmodel!(v, x, θ)

        end

    end
    
    # Levenberg-Marquardt algorithm
    # -----------------------------
    
    Id = Matrix(1.0I, n, n)
    
    # Status = 1 means success
    status = 1

    # Parameters
    
    λ_up      = 2.0
    λ_down    = 2.0
    λ         = 1.0
    maxoutind = min(p, 5)

    # Allocation
    θnew = Vector{Float64}(undef, n)
    d = Vector{Float64}(undef, n)
    y = Vector{Float64}(undef, n)
    G = Array{Float64, 2}(undef, n, n)
    grad_lovo = Vector{Float64}(undef, n)
    
    # Initial steps
    
    ind_lovo, best_val_lovo = LovoFun(θ)

    ResFun!(θ, ind_lovo, val_res, jac_res)

    BLAS.gemv!('T', 1.0, jac_res, val_res, 0.0, grad_lovo)
 
    ngrad_lovo = norm(grad_lovo, 2)
    
    safecount = 1

    # Main loop
    
    while (ngrad_lovo >= ε) && (safecount < MAXITER)

        with_logger(lm_logger.x) do

            @info("Iteration $(safecount)")
            @info("  Current value:   $(best_val_lovo)")
            @info("  ||grad_lovo||_2: $(ngrad_lovo)")
            @info("  Current iterate: $(θ)")
            @info("  Best indices (first $(maxoutind)): $(ind_lovo[1:maxoutind])")
            @info("  lambda: $(λ)")

        end

        G .= Id

        BLAS.gemm!('T', 'N', 1.0, jac_res, jac_res, λ, G)
        
        F = qr!(G)
        #F = cholesky!(G, Val(true))
        
        ad = try

            ldiv!(d, F, grad_lovo)

            d .*= -1.0
            
        catch
            "error"
        end
        if ad == "error" #restarting if lapack fails
            with_logger(lm_logger.x) do
            
                @warn "Failed to solve the linear system. Will try new point."
                d .= - 1.0 .* grad_lovo 
                θ .= rand(n)
                
            end
        else 
            d .= ad
        end

        θnew .= θ .+ d
        
        ind_lovo, val_lovo = LovoFun(θnew)
        
        if  val_lovo <= best_val_lovo

            θ .= θnew
            
            best_val_lovo = val_lovo
            
            λ = λ / λ_down
            
            ResFun!(θ, ind_lovo, val_res, jac_res)

            BLAS.gemv!('T', 1.0, jac_res, val_res, 0.0, grad_lovo)

            ngrad_lovo = norm(grad_lovo, 2)

            with_logger(lm_logger.x) do
                
                @info("  Better function value found, lambda changed to $(λ).")

            end
            
        else

            λ = λ * λ_up

            with_logger(lm_logger.x) do

                @info("  No improvement, lambda changed to $(λ).")

            end
            
        end

        safecount += 1
        
    end
    
    if safecount == MAXITER
        with_logger(lm_logger.x) do
            
            @info("No solution was found in $(safecount) iterations.")

        end
        status = 0
    end

    # TODO: Create a test for this case
    if isnan(ngrad_lovo)
        with_logger(lm_logger.x) do

            @info("Incorrect value for gradient norm $(ngrad_lovo).")

        end
        status = 0
    end

    outliers = [1:npun;]
    setdiff!(outliers, ind_lovo)

    with_logger(lm_logger.x) do

        @info("""

        Final iteration (STATUS=$(status))
          Solution found:       $(θ)
          ||grad_lovo||_2:      $(ngrad_lovo)
          Function value:       $(best_val_lovo)
          Number of iterations: $(safecount)
          Outliers:             $(outliers)
    
        """)

    end
    
    return RAFFOutput(status, θ, safecount, p, best_val_lovo, nf, nj, outliers)

end

function lmlovo(model::Function, θ::Vector{Float64}, data::Array{Float64,2},
                n::Int, p::Int; kwargs...)

    # Define closures for derivative and initializations

    # 'x' is considered as global parameter for this function
    model_cl(θ) = model(x, θ)
    
    grad_model!(g, x_, θ) = begin

        global x = x_
        
        ForwardDiff.gradient!(g, model_cl, θ)

    end

    return lmlovo(model, grad_model!, θ, data, n, p; kwargs...)
    
end

lmlovo(model::Function, gmodel!::Function, data::Array{Float64,2}, n::Int, p::Int; kwargs...) =
    lmlovo(model, gmodel!, zeros(Float64, n), data, n, p; kwargs...)

lmlovo(model::Function, data::Array{Float64,2}, n::Int, p::Int; kwargs...) =
    lmlovo(model, zeros(Float64, n), data, n, p; kwargs...)

"""

    raff(model::Function, data::Array{Float64, 2}, n::Int; kwargs...)

    raff(model::Function, gmodel!::Function, data::Array{Float64, 2},
        n::Int; MAXMS::Int=1, SEEDMS::Int=123456789,
        initguess::Vector{Float64}=zeros(Float64, n),
        ε::Float64=1.0e-4, noutliers::Int=-1, ftrusted::Union{Float64,
        Tuple{Float64, Float64}}=0.5)

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
  - `ε`: gradient stopping criteria to `lmlovo`
  - `noutliers`: integer describing the maximum expected number of
    outliers. The default is *half*. *Deprecated*.
  - `ftrusted`: float describing the minimum expected percentage of
    trusted points. The default is *half* (0.5). Can also be a
    Tuple of the form `(fmin, fmax)` percentages of trusted points.

Returns a [`RAFFOutput`](@ref) object with the best parameter found.

"""
function raff(model::Function, gmodel!::Function,
    data::Array{Float64, 2}, n::Int; MAXMS::Int=1,
    SEEDMS::Int=123456789, initguess::Vector{Float64}=zeros(Float64,
    n), ε::Float64=1.0e-4, noutliers::Int=-1,
    ftrusted::Union{Float64, Tuple{Float64, Float64}}=0.5)

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
                
                @debug("Running lmlovo for p = $(i). Repetition $(j).")
                
            end
            
            # Starting point
            θ = randn(seedMS, Float64, n)
            θ .= θ .+ initguess
        
            # Call function and store results
            sols[ind] = lmlovo(model, gmodel!, θ, data, n, i; ε=ε, MAXITER=400)

            # Update the best point and functional value
            (sols[ind].status == 1) && (sols[ind].f < vbest.f) && (vbest = sols[ind])
            
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
        ε::Float64=1.0e-4, noutliers::Int=-1, ftrusted::Union{Float64,
        Tuple{Float64, Float64}}=0.5)

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
  - `ε`: stopping tolerance
  - `noutliers`: integer describing the maximum expected number of
    outliers. The default is *half*. *Deprecated*.
  - `ftrusted`: float describing the minimum expected percentage of
    trusted points. The default is *half* (0.5). Can also be a
    Tuple of the form `(fmin, fmax)` percentages of trusted points.

Returns a [`RAFFOutput`](@ref) object containing the solution.

"""
function praff(model::Function, gmodel!::Function,
               data::Array{Float64, 2}, n::Int; MAXMS::Int=1,
               SEEDMS::Int=123456789, batches::Int=1,
               initguess::Vector{Float64}=zeros(Float64, n),
               ε::Float64=1.0e-4, noutliers::Int=-1,
               ftrusted::Union{Float64, Tuple{Float64, Float64}}=0.5)

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
                           plimsup, MAXMS, seedMS, initguess)
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
