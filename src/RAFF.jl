__precompile__(false)

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

# Set RAFF logger
raff_logger = ConsoleLogger(stdout, Logging.Error)

lm_logger = ConsoleLogger(stdout, Logging.Error)

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

    t11 t12 ... t1N y1
    t21 t22 ... t2N y2
    :

where `N` is the dimension of the argument of the model
(i.e. dimension of `t`).

If 'θ' is provided, then it is used as the starting point.

The signature of function `model` should be given by

    model(x::Union{Vector{Float64}, SubArray}, θ::Vector{Float64})

where `x` are the variables and `θ` is a `n`-dimensional vector of
parameters. If the gradient of the model `gmodel!`

    gmodel!(g::Vector{Float64}, x::Union{Vector{Float64},
            θ::Vector{Float64}, SubArray} )

is not provided, then the function ForwardDiff.gradient! is called to
compute it.  **Note** that this choice has an impact in the
computational performance of the algorithm. In addition, if
`ForwardDiff.jl` is being used, then one **MUST** remove the signature
of vector `θ` from the model.

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

    with_logger(lm_logger) do
        
        @debug("Size of data matrix ", size(data))

    end

    (p == 0) && return RAFFOutput(1, θ, 0, p, 0.0, [1:npun;])
    
    # Main function - the LOVO function
    LovoFun = let

        npun_::Int = npun
        
        ind::Vector{Int} = Vector{Int}(undef, npun_)
        
        F::Vector{Float64} = Vector{Float64}(undef, npun_)
        
        p_::Int = p
        
        # Return a ordered set index and lovo value
        (θ) -> begin
            
            @views for i = 1:npun_
                F[i] = (model(data[i,1:(end - 1)]) - data[i, end], θ)^2
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

        with_logger(lm_logger) do

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
            with_logger(lm_logger) do
            
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

            with_logger(lm_logger) do
                
                @info("  Better function value found, lambda changed to $(λ).")

            end
            
        else

            λ = λ * λ_up

            with_logger(lm_logger) do

                @info("  No improvement, lambda changed to $(λ).")

            end
            
        end

        safecount += 1
        
    end
    
    if safecount == MAXITER
        with_logger(lm_logger) do
            
            @info("No solution was found in $(safecount) iterations.")

        end
        status = 0
    end

    # TODO: Create a test for this case
    if isnan(ngrad_lovo)
        with_logger(lm_logger) do

            @info("Incorrect value for gradient norm $(ngrad_lovo).")

        end
        status = 0
    end

    outliers = [1:npun;]
    setdiff!(outliers, ind_lovo)

    with_logger(lm_logger) do

        @info("""

        Final iteration (STATUS=$(status))
          Solution found:       $(θ)
          ||grad_lovo||_2:      $(ngrad_lovo)
          Function value:       $(best_val_lovo)
          Number of iterations: $(safecount)
          Outliers:             $(outliers)
    
        """)

    end
    
    return RAFFOutput(status, θ, safecount, p, best_val_lovo, outliers)

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

    return lmlovo(model, grad_model, θ, data, n, p; kwargs...)
    
end

lmlovo(model::Function, gmodel!::Function, data::Array{Float64,2}, n::Int, p::Int; kwargs...) =
    lmlovo(model, gmodel!, zeros(Float64, n), data, n, p; kwargs...)

lmlovo(model::Function, data::Array{Float64,2}, n::Int, p::Int; kwargs...) =
    lmlovo(model, zeros(Float64, n), data, n, p; kwargs...)

"""
    raff(model::Function, data::Array{Float64, 2}, n::Int; MAXMS::Int=1,
         SEEDMS::Int=123456789, initguess=zeros(Float64, n))

    raff(model::Function, gmodel!::Function, data::Array{Float64, 2}, n::Int;
         [MAXMS::Int=1, SEEDMS::Int=123456789, initguess=zeros(Float64, n),
          kwargs...])

Robust Algebric Fitting Function (RAFF) algorithm. This function uses
a voting system to automatically find the number of trusted data
points to fit the `model`.

  - `model`: function to fit data. Its signature should be given by

        model(x, t)

    where `x` is a `n`-dimensional vector of parameters and `t` is the
    multidimensional argument

  - `gmodel!`: gradient of the model function. Its signature should be
    given by

        gmodel!(x, t, g)

    where `x` is a `n`-dimensional vector of parameters, `t` is the
    multidimensional argument and the gradient is written in `g`.

  - `data`: data to be fit. This matrix should be in the form

        t11 t12 ... t1N y1
        t21 t22 ... t2N y2
        :

    where `N` is the dimension of the argument of the model
    (i.e. dimension of `t`).

  - `n`: dimension of the parameter vector in the model function

The optional arguments are

  - `MAXMS`: number of multistart points to be used
  - `SEEDMS`: integer seed for random multistart points
  - `initialguess`: a good guess for the starting point and for
    generating random points in the multistart strategy
  - `ε`: gradient stopping criteria to `lmlovo`
  - `noutliers`: integer describing the maximum expected number of
    outliers. The default is *half*.

Returns a [`RAFFOutput`](@ref) object with the best parameter found.

"""
function raff(model::Function, gmodel!::Function,
              data::Array{Float64, 2}, n::Int; MAXMS::Int=1,
              SEEDMS::Int=123456789, initguess=zeros(Float64, n), ε=1.0e-4,
              noutliers::Int=-1)

    # Initialize random generator
    seedMS = MersenneTwister(SEEDMS)

    plimsup = length(data[:, 1])
    pliminf = (noutliers >= 0) ? plimsup - noutliers : Int(round(length(data[:, 1]) / 2.0))
    lv = plimsup - pliminf + 1
    
    sols = Vector{RAFFOutput}(undef, lv)

    for i = pliminf:plimsup

        vbest = RAFFOutput(0, initguess, -1, i, Inf, [])
        
        ind = i - pliminf + 1
        
        for j = 1:MAXMS

            with_logger(raff_logger) do
                
                @debug("Running lmlovo for p = $(i). Repetition $(j).")
                
            end
            
            # Starting point
            x = randn(seedMS, Float64, n)
            # x .= x .+ vbest.solution
        
            # Call function and store results
            sols[ind] = lmlovo(model, gmodel!, x, data, n, i; ε=ε, MAXITER=400)

            # Update the best point and functional value
            (sols[ind].status == 1) && (sols[ind].f < vbest.f) && (vbest = sols[ind])
            
        end

        sols[ind] = vbest

    end

    # Remove possible stationary points, i.e., points with lower
    # values for 'p' and higher 'f'.

    with_logger(raff_logger) do
        
        eliminate_local_min!(model, data, sols)

    end

    # Voting strategy
    
    dvector = zeros(Int(lv * (lv - 1) / 2))
    dmatrix = zeros(lv, lv)
    pos = 0
    n_conv = 0

    for j = 1:lv

        # Count how many have successfully converged
        (sols[j].status == 1) && (n_conv += 1)
        
        for i = j + 1:lv

            dmatrix[i, j] = Inf

            if sols[i].status == 1 && sols[j].status == 1

                dmatrix[i, j] = norm(sols[i].solution - sols[j].solution)

                pos += 1

                dvector[pos] = dmatrix[i, j]

            end

        end

    end

    threshold = Inf

    if pos > 0
        
        dvv = @view dvector[1:pos]
        
        threshold = minimum(dvv) + mean(dvv) / (1.0 + sqrt(plimsup))

    elseif n_conv == 0

        with_logger(raff_logger) do
            
            @warn("No convergence for any 'p'. Returning largest.")

        end

    end
    
    votsis = zeros(lv)
    with_logger(raff_logger) do;  @debug("Threshold: $(threshold)"); end

    # Actual votation
    
    for j = 1:lv
        # Count +1 if converged
        (sols[j].status == 1) && (votsis[j] += 1)
        # Check other distances
        for i = j + 1:lv
            if dmatrix[i, j] <=  threshold
                votsis[j] += 1
                votsis[i] += 1
            end
        end
    end

    with_logger(raff_logger) do

        @debug("Voting vector:", votsis)
        @debug("Distance matrix:", dmatrix)

    end
    
    mainind = findlast(x->x == maximum(votsis), votsis)
    
    return sols[mainind]
    
end

function raff(model::Function, data::Array{Float64, 2}, n::Int; kwargs...)

    # Define closures for derivative and initializations

    # 't' is considered as global parameter for this function
    model_cl(x) = model(x, t)
    
    grad_model(x, t_, g) = begin

        global t = t_
        
        return ForwardDiff.gradient!(g, model_cl, x)

    end

    return raff(model, grad_model, data, n; kwargs...)

end

"""
    praff(model::Function, data::Array{Float64, 2}, n::Int; MAXMS::Int=1,
          SEEDMS::Int=123456789, batches::Int=1, initguess=zeros(Float64, n),
          ε=1.0e-4)

    praff(model::Function, gmodel!::Function, data::Array{Float64, 2}, n::Int;
          MAXMS::Int=1, SEEDMS::Int=123456789, batches::Int=1,
          initguess=zeros(Float64, n), ε::Float64=1.0e-4)

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
    outliers. The default is *half*.

Returns a [`RAFFOutput`](@ref) object containing the solution.

"""
function praff(model::Function, gmodel!::Function,
               data::Array{Float64, 2}, n::Int; MAXMS::Int=1,
               SEEDMS::Int=123456789, batches::Int=1, initguess=zeros(Float64, n),
               ε::Float64=1.0e-4, noutliers::Int=-1)

    # Initialize random generator
    seedMS = MersenneTwister(SEEDMS)
    
    plimsup = length(data[:, 1])
    pliminf = (noutliers >= 0) ? plimsup - noutliers : Int(round(length(data[:, 1]) / 2.0))
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
    with_logger(raff_logger) do

        @debug("Workers", curr_workers)

    end
    
    for (i, t) in enumerate(curr_workers)

        futures[i] = @spawnat(t, with_logger( ()-> try

            @debug("Creating worker $(t).")
                                              
            consume_tqueue(bqueue, tqueue, squeue,
                           model, gmodel!, data, n, pliminf,
                           plimsup, MAXMS, seedMS)
            catch e
                                              
               @error("Unable to start worker $(t).", e)

            end, raff_logger
        ))

        with_logger(raff_logger) do

            @debug("Created worker $(t).")

        end
    end

    # Check asynchronously if there is at least one live worker
    @async with_logger(
        () -> check_and_close(bqueue, tqueue, squeue, futures),
        raff_logger)

    # Populate the task queue with jobs
    for p = pliminf:batches:plimsup

        try
            
            put!(tqueue, p:min(plimsup, p + batches - 1))

        catch e

            with_logger(raff_logger) do

                @warn("Tasks queue prematurely closed while inserting tasks. Will exit.")

            end

            break

        end

        with_logger(raff_logger) do

            @debug("Added problem $(p) to tasks queue.")

        end
                
    end

    # The task queue can be closed, since all the problems have been
    # read, due to the size 0 of this channel
    close(tqueue)

    with_logger(raff_logger) do
        
        @debug("Waiting for workers to finish.")

    end

    # Create a vector of solutions to store the results from workers
    sols = Vector{RAFFOutput}(undef, lv)

    for i in 1:lv

        try

            rout = take!(squeue)

            sols[rout.p - pliminf + 1] = rout

            with_logger(raff_logger) do

                @debug("Stored solution for p=$(rout.p).")

            end

        catch e

            with_logger(raff_logger) do

                @error("Error when retrieving solutions.", e)

            end
            
        end
        
    end
    
    close(bqueue)

    close(squeue)

    # Voting strategy

    dvector = zeros(Int(lv * (lv - 1) / 2))
    dmatrix = zeros(lv, lv)
    pos = 0
    n_conv = 0

    for j = 1:lv

        # Count how many have successfully converged
        (sols[j].status == 1) && (n_conv += 1)

        for i = j + 1:lv

            dmatrix[i, j] = Inf

            if sols[i].status == 1 && sols[j].status == 1

                dmatrix[i, j] = norm(sols[i].solution - sols[j].solution)

                pos += 1

                dvector[pos] = dmatrix[i, j]

            end

        end

    end

    threshold = Inf

    if pos > 0

        dvv = @view dvector[1:pos]

        threshold = minimum(dvv) + mean(dvv) / (1.0 + sqrt(plimsup))

    elseif n_conv == 0

        with_logger(raff_logger) do
            
            @warn("No convergence for any 'p'. Returning largest.")

        end

    end
    
    votsis = zeros(lv)

    with_logger(raff_logger) do;  @debug("Threshold: $(threshold)"); end

    # Actual votation
    
    for j = 1:lv
        # Count +1 if converged
        (sols[j].status == 1) && (votsis[j] += 1)
        # Check other distances
        for i = j + 1:lv
            if dmatrix[i, j] <=  threshold
                votsis[j] += 1
                votsis[i] += 1
            end
        end
    end

    with_logger(raff_logger) do

        @debug("Voting vector:", votsis)
        @debug("Distance matrix:", dmatrix)

    end
    
    mainind = findlast(x->x == maximum(votsis), votsis)
    
    return sols[mainind]

end

function praff(model::Function, data::Array{Float64, 2}, n::Int; kwargs...)

    # Define closures for derivative and initializations

    # 't' is considered as global parameter for this function
    model_cl(x) = model(x, t)
    
    grad_model(x, t_, g) = begin

        global t = t_
        
        return ForwardDiff.gradient!(g, model_cl, x)

    end

    return praff(model, grad_model, data, n; kwargs...)

end


end
