module RAFF

__precompile__(false)

# Dependencies
using Distributed
using ForwardDiff
using LinearAlgebra
using Statistics
using Printf
using Random
using SharedArrays

export lmlovo, raff, praff

include("raffoutput.jl")

include("utils.jl")

include("dutils.jl")

include("generator.jl")

"""
    lmlovo(model::Function [, x::Vector{Float64} = zeros(n)], data::Array{Float64, 2},
           n::Int, p::Int [; kwargs...])

    lmlovo(model::Function, gmodel!::Function [, x::Vector{Float64} = zeros(n)],
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

If 'x' is provided, the it is used as the starting point.

The signature of function `model` should be given by

    model(x::Vector{Float64}, t::Union{Vector{Float64}, SubArray})

where `x` is a `n`-dimensional vector of parameters and `t` is the
argument. If the gradient of the model `gmodel!`

    gmodel!(x::Vector{Float64}, t::Union{Vector{Float64}, SubArray},
            g::Vector{Float64})

is not provided, then the function ForwardDiff.gradient! is called to
compute it.  **Note** that this choice has an impact in the
computational performance of the algorithm. In addition, if
ForwardDiff is being used, then one **MUST** remove the signature of
vector `x` from the model.

The optional arguments are

  - `MAXITER`: maximum number of iterations
  - `ε`: tolerance for the gradient of the function

Returns a RAFFOutput object.

"""
function lmlovo(model::Function, gmodel!::Function, x::Vector{Float64},
                data::Array{Float64,2}, n::Int, p::Int;
                MAXITER::Int=200, ε::Float64=10.0^-4)

    @assert(n > 0, "Dimension should be positive.")
    @assert(p >= 0, "Trusted points should be nonnegative.")
    
    npun, = size(data)

    @debug("Size of data matrix ", size(data))
    
    # Main function - the LOVO function
    LovoFun = let

        npun_::Int = npun
        
        ind::Vector{Int} = Vector{Int}(undef, npun_)
        
        F::Vector{Float64} = Vector{Float64}(undef, npun_)
        
        p_::Int = p
        
        # Return a ordered set index and lovo value
        x -> begin
            
            @views for i = 1:npun_
                F[i] = (model(x, data[i,1:(end - 1)]) - data[i, end])^2
            end
            
            indF, orderedF = SortFun!(F, ind, p_)
            
            return indF, sum(orderedF)
        end
        
    end
    
    # Residue and Jacobian of residue
    val_res::Vector{Float64} = Vector{Float64}(undef, p)
    
    jac_res::Array{Float64, 2} = Array{Float64}(undef, p, n)
    
    # This function returns the residue and Jacobian of residue
    ResFun!(x::Vector{Float64}, ind, r::Vector{Float64},
            rJ::Array{Float64, 2}) = begin
                
        for (k, i) in enumerate(ind)
            
            t = data[i, 1:(end - 1)]
            
            r[k] = model(x, t) - data[i, end]
            
            v = @view(rJ[k, :])
            
            gmodel!(x, t, v)

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
    xnew = Vector{Float64}(undef, n)
    d = Vector{Float64}(undef, n)
    y = Vector{Float64}(undef, n)
    G = Array{Float64, 2}(undef, n, n)
    grad_lovo = Vector{Float64}(undef, n)
    
    # Initial steps
    
    ind_lovo, best_val_lovo = LovoFun(x)

    ResFun!(x, ind_lovo, val_res, jac_res)

    BLAS.gemv!('T', 1.0, jac_res, val_res, 0.0, grad_lovo)
 
    ngrad_lovo = norm(grad_lovo, 2)
    
    safecount = 1

    # Main loop
    
    while (ngrad_lovo >= ε) && (safecount < MAXITER)

        @info("Iteration $(safecount)")
        @info("  Current value:   $(best_val_lovo)")
        @info("  ||grad_lovo||_2: $(ngrad_lovo)")
        @info("  Current iterate: $(x)")
        @info("  Best indices (first $(maxoutind)): $(ind_lovo[1:maxoutind])")
        @info("  lambda: $(λ)")

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
            @warn "Failed to solve the linear system. Will try new point."
            d .= - 1.0 .* grad_lovo 
            x .= rand(n)
        else 
            d .= ad
        end

        xnew .= x .+ d
        
        ind_lovo, val_lovo = LovoFun(xnew)
        
        if  val_lovo <= best_val_lovo

            x .= xnew
            
            best_val_lovo = val_lovo
            
            λ = λ / λ_down
            
            ResFun!(x, ind_lovo, val_res, jac_res)

            BLAS.gemv!('T', 1.0, jac_res, val_res, 0.0, grad_lovo)

            ngrad_lovo = norm(grad_lovo, 2)
            
            @info("  Better function value found, lambda changed to $(λ).")
            
        else

            λ = λ * λ_up
            
            @info("  No improvement, lambda changed to $(λ).")
            
        end

        safecount += 1
        
    end
    
    if safecount == MAXITER
        @info("No solution was found in $(safecount) iterations.")
        status = 0
    end

    # TODO: Create a test for this case
    if isnan(ngrad_lovo)
        @info("Incorrect value for gradient norm $(ngrad_lovo).")
        status = 0
    end

    outliers = [1:npun;]
    setdiff!(outliers, ind_lovo)
    
    @info("""

    Final iteration (STATUS=$(status))
      Solution found:       $(x)
      ||grad_lovo||_2:      $(ngrad_lovo)
      Function value:       $(best_val_lovo)
      Number of iterations: $(safecount)
      Outliers:             $(outliers)
    
    """)
    
    return RAFFOutput(status, x, safecount, p, best_val_lovo, outliers)

end

function lmlovo(model::Function, x::Vector{Float64}, data::Array{Float64,2},
                n::Int, p::Int; kwargs...)

    # Define closures for derivative and initializations

    # 't' is considered as global parameter for this function
    model_cl(x) = model(x, t)
    
    grad_model(x, t_, g) = begin

        global t = t_
        
        ForwardDiff.gradient!(g, model_cl, x)

    end

    return lmlovo(model, grad_model, x, data, n, p; kwargs...)
    
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

Returns a RAFFOutput object with the best parameter found.

"""
function raff(model::Function, gmodel!::Function,
              data::Array{Float64, 2}, n::Int; MAXMS::Int=1,
              SEEDMS::Int=123456789, initguess=zeros(Float64, n), ε=1.0e-4)

    # Initialize random generator
    seedMS = MersenneTwister(SEEDMS)

    pliminf = Int(round(length(data[:, 1]) / 2.0))
    plimsup = length(data[:, 1])
    lv = plimsup - pliminf + 1
    
    sols = Vector{RAFFOutput}(undef, lv)

    for i = pliminf:plimsup

        vbest = RAFFOutput(0, initguess, -1, i, Inf, [])
        
        ind = i - pliminf + 1
        
        for j = 1:MAXMS
            
            @debug("Running lmlovo for p = $(i). Repetition $(j).")
        
            # Starting point
            x = randn(seedMS, Float64, n)
            x .= x .+ vbest.solution
        
            # Call function and store results
            sols[ind] = lmlovo(model, gmodel!, x, data, n, i; ε=ε)

            # Update the best point and functional value
            (sols[ind].status == 1) && (sols[ind].f < vbest.f) && (vbest = sols[ind])
            
        end

        sols[ind] = vbest

    end

    # Remove possible stationary points, i.e., points with lower
    # values for 'p' and higher 'f'.

    eliminate_local_min!(model, data, sols)

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
        
        @warn("No convergence for any 'p'. Returning largest.")

    end
    
    votsis = zeros(lv)
    @debug("Threshold: $(threshold)")

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
    
    @debug("Voting vector:", votsis)
    @debug("Distance matrix:", dmatrix)
    
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

    praff(model::Function, data::Array{Float64, 2}, n::Int;
          MAXMS::Int=1, SEEDMS::Int=123456789, batches::Int=1 )

    praff(model::Function, gmodel!::Function,
          data::Array{Float64, 2}, n::Int; MAXMS::Int=1,
          SEEDMS::Int=123456789, batches::Int=1 )

Multicore shared memory version of RAFF. See the description of the
[raff](@ref) function for the main (non-optional) arguments.

This function uses all available **local** workers to run RAFF
algorithm. Note that this function does not use *Tasks*, so all the
parallelism is based on the `Distributed` package.

The optional arguments are

  - `MAXMS`: number of multistart points to be used
  - `SEEDMS`: integer seed for random multistart points
  - `batches`: size of batches to be send to each worker

Returns a tuple `x`, `f`, `p` where

  - `x` is the solution
  - `f` is the value of the error at the solution
  - `p` is the number of trusted points found by the voting system.

"""
function praff(model::Function, gmodel!::Function,
               data::Array{Float64, 2}, n::Int; MAXMS::Int=1,
               SEEDMS::Int=123456789, batches::Int=1)

    # Initialize random generator
    seedMS = MersenneTwister(SEEDMS)
    
    pliminf = Int(round(length(data[:, 1]) / 2.0))
    plimsup = length(data[:, 1])

    nInfo = plimsup - pliminf + 1
    
    v = SharedArray{Float64, 2}(n, nInfo)
    vf = SharedArray{Float64, 1}(nInfo)
    vs = SharedArray{Int, 1}(nInfo)
    bestx = SharedArray{Float64, 1}(n)

    vs .= 0
    
    # Create a RemoteChannel to receive solutions
    bqueue = RemoteChannel(() -> Channel{Vector{Float64}}(div(nInfo, 2)))
    # Create another channel to assign tasks
    tqueue = RemoteChannel(() -> Channel{UnitRange{Int}}(0))

    # This command selects only nodes which are local to myid()
    local_workers = intersect(workers(), procs(myid()))
    futures = Vector{Future}(undef, length(local_workers))
    
    # Start updater Task
    @async update_best(bqueue, bestx)

    # Start workers Tasks (CPU intensive)
    for (i, t) in enumerate(local_workers)

        futures[i] = @spawnat(t, consume_tqueue( bqueue, tqueue,
            bestx, v, vs, vf, model, gmodel!, data, n, pliminf,
            plimsup, MAXMS, seedMS ))
        
    end

    # Check asynchronously if there is at least one live worker
    @async check_and_close(bqueue, tqueue, futures)

    # Populate the task queue with jobs
    for p = pliminf:batches:plimsup

        try
            
            put!(tqueue, p:min(plimsup, p + batches))

        catch e

            @warn("Tasks queue prematurely closed while inserting tasks. Will exit.")

            break

        end

        @debug("Added problem $(p) to tasks queue.")
                
    end

    # The task queue can be closed, since all the problems have been
    # read, due to the size 0 of this channel
    close(tqueue)

    @debug("Waiting for workers to finish.")

    for f in futures

        try

            wait(f)

        catch e

            @error("Error in consumer for worker $(f.where)", e)

        end
        
    end
    
    close(bqueue)
    
    votsis = zeros(nInfo)
    
    for i = 1:nInfo
        for j = 1:nInfo
            if vs[i] == 1 && vs[j] == 1 && norm(v[:, i] - v[:, j]) < 10.0^(-3)
                votsis[i] += 1
            end
        end
    end
    
    mainind = findlast(x->x == maximum(votsis), votsis)

    @info("""Solution from PRAFF:
    x= $(v[:, mainind])
    f= $(vf[mainind])
    p= $(pliminf + mainind - 1)
    """)
    
    return v[:, mainind], vf[mainind], pliminf + mainind - 1
    
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
