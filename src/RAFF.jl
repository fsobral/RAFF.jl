module RAFF

__precompile__(false)

# Dependencies
using Distributed
using ForwardDiff
using LinearAlgebra
using Printf
using Random
using SharedArrays

export LMlovo, raff, praff, generateTestProblems

include("utils.jl")

"""
    LMlovo(model::Function [, x::Vector{Float64} = zeros(n)], data::Array{Float64, 2},
           n::Int, p::Int [; kwargs...])

    LMlovo(model::Function, gmodel!::Function [, x::Vector{Float64} = zeros(n)],
           data::Array{Float64,2}, n::Int, p::Int [; MAXITER::Int])

Fit the `n`-parameter model `model` to the data given by matrix
`data`. The strategy is based on the LOVO function, which means that
only `p` (0 < `p` <= `n`) points are trusted. The Levenberg-Marquardt
algorithm is implemented in this version.

If 'x' is provided, the it is used as the starting point.

The signature of function `model` should be given by

    model(x::Vector{Float64}, t::Float64)

where `x` is a `n`-dimensional vector of parameters and `t` is the
argument. If the gradient of the model `gmodel!`

    gmodel!(x::Vector{Float64}, t::Float64, g::Vector{Float64})

is not provided, then the function ForwardDiff.gradient! is called to
compute it. **Note** that this choice has an impact in the
computational performance of the algorithm.

Returns a tuple `s`, `x`, `iter`, `p`, where

  - `s`: is 1 if converged and 0 if not
  - `x`: vector with the parameters of the model
  - `iter`: number of iterations up to convergence
  - `p`: number of trusted points
  - `f`: the residual value

"""
function LMlovo(model::Function, gmodel!::Function, x::Vector{Float64},
                data::Array{Float64,2}, n::Int, p::Int;
                MAXITER::Int=200)

    @assert(n > 0, "Dimension should be positive.")
    @assert(p >= 0, "Trusted points should be nonnegative.")
    
    npun = length(data[:, 1])
    
    # Main function - the LOVO function
    LovoFun = let

        npun_::Int = npun
        
        ind::Vector{Int} = Vector{Int}(undef, npun_)

        F::Vector{Float64} = Vector{Float64}(undef, npun_)

        p_::Int = p

        # Return a ordered set index and lovo value
        x -> begin
        
            for i = 1:npun_
                F[i] = (model(x, data[i,1]) - data[i, 2])^2
            end
            
            indF, orderedF = SortFun!(F, ind, p_)
            
            return indF, sum(orderedF)
        end

    end

    # Residue and Jacobian of residue
    val_res::Vector{Float64} = Vector{Float64}(undef, p)

    jac_res::Array{Float64, 2} = Array{Float64}(undef, p, n)
    
    # This function returns the residue and Jacobian of residue 
    ResFun!(x::Vector{Float64}, ind,
            r::Vector{Float64}, rJ::Array{Float64, 2}) = begin

        k = 1
                
        for i in ind
            t = data[i, 1]
            r[k] = model(x, t) - data[i, 2]

            v = @view(rJ[k, :])
            gmodel!(x, t, v)

            k = k + 1
        end

    end
    
    # Levenberg-Marquardt algorithm
    # -----------------------------
    
    Id = Matrix(1.0I, n, n)

    # Status = 1 means success
    status = 1

    # Parameters
    ε         = 10.0^(-4)
    λ_up      = 2.0
    λ_down    = 2.0
    λ         = 1.0
    MAXOUTIND = 5

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
        @info("  Best indices (first $(MAXOUTIND)): $(ind_lovo[1:MAXOUTIND])")
        @info("  lambda: $(λ)")
        
        G .= Id
        
        BLAS.gemm!('T', 'N', 1.0, jac_res, jac_res, λ, G)
        
        F = qr!(G)
        
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

    @info("""

    Final iteration (STATUS=$(status))
      Solution found:       $(x)
      ||grad_lovo||_2:      $(ngrad_lovo)
      Function value:       $(best_val_lovo)
      Number of iterations: $(safecount)

    """)
    
    return status, x, safecount, p, best_val_lovo

end

function LMlovo(model::Function, x::Vector{Float64}, data::Array{Float64,2},
                n::Int, p::Int; kwargs...)

    # Define closures for derivative and initializations

    # 't' is considered as global parameter for this function
    model_cl(x) = model(x, t)
    
    grad_model(x, t_, g) = begin

        global t = t_
        
        ForwardDiff.gradient!(g, model_cl, x)

    end

    return LMlovo(model, grad_model, x, data, n, p; kwargs...)
    
end

LMlovo(model::Function, gmodel!::Function, data::Array{Float64,2}, n::Int, p::Int; kwargs...) =
    LMlovo(model, gmodel!, zeros(Float64, n), data, n, p; kwargs...)

LMlovo(model::Function, data::Array{Float64,2}, n::Int, p::Int; kwargs...) =
    LMlovo(model, zeros(Float64, n), data, n, p; kwargs...)

"""
    raff(model::Function, data::Array{Float64, 2}, n::Int)

    raff(model::Function, gmodel!::Function, data::Array{Float64, 2}, n::Int)

Robust Algebric Fitting Function (RAFF) algorithm. This function uses
a voting system to automatically find the number of trusted data
points to fit the `model`.

  - `model`: function to fit data. Its signature should be given by

      model(x, t)

    where `x` is a `n`-dimensional vector of parameters and `t` is the
    argument
  - `gmodel!`: gradient of the model function. Its signature should be
    given by

      gmodel!(x, t, g)

    where `x` is a `n`-dimensional vector of parameters and `t` is the
    argument and the gradient is written in `g`.
  - `data`: data to be fit
  - `n`: dimension of the parameter vector in the model function

"""
function raff(model::Function, gmodel!::Function,
              data::Array{Float64, 2}, n::Int)

    pliminf = Int(round(length(data[:, 1]) / 2.0))
    plimsup = length(data[:, 1])
    
    v = Array{Any,1}(undef, plimsup - pliminf + 1)

    for i = pliminf:plimsup

        @debug("Running LMlovo for p = $(i).")
        
        # Starting point
        x = zeros(Float64, n)
        
        # Call function and store results
        v[i - pliminf + 1] = LMlovo(model, gmodel!, x, data, n, i)
        
    end
    
    lv = length(v)
    votsis = zeros(lv)
    for i = 1:lv
        for j = 1:lv
            if norm(v[i][2] - v[j][2]) < 10.0^(-3)
                votsis[i] += 1
            end
        end
    end
    
    mainind = findlast(x->x == maximum(votsis), votsis)
    
    return v[mainind]
    
end

function raff(model::Function, data::Array{Float64, 2}, n::Int)

    # Define closures for derivative and initializations

    # 't' is considered as global parameter for this function
    model_cl(x) = model(x, t)
    
    grad_model(x, t_, g) = begin

        global t = t_
        
        return ForwardDiff.gradient!(g, model_cl, x)

    end

    return raff(model, grad_model, data, n)

end

"""

    praff(model::Function, gmodel!::Function,
            data::Array{Float64, 2}, n::Int)

Parallel and shared memory version of RAFF. See the description of the
[raff](@ref) function.

This function uses all available local workers to run the RAFF
algorithm.

"""
function praff(model::Function, gmodel!::Function,
               data::Array{Float64, 2}, n::Int; MAXMS::Int=1,
               SEEDMS::Int=1234)

    # Initializes random generator
    seed = MersenneTwister(SEEDMS)
    
    pliminf = Int(round(length(data[:, 1]) / 2.0))
    plimsup = length(data[:, 1])

    nInfo = plimsup - pliminf + 1
    
    v = SharedArray{Float64, 2}(n, nInfo)
    vf = SharedArray{Float64, 1}(nInfo)

    f = @sync @distributed for i = pliminf:plimsup

        # Starting point
        bestx = zeros(Float64, n)
        bestf = Inf

        # Multi-start strategy
        for j = 1:MAXMS

            # New random starting point
            x  = randn(seed, n)
            x .= 10.0 .* x .+ bestx
        
            # Call function and store results
            s, x, iter, p, f = LMlovo(model, gmodel!, x, data, n, i)

            if f < bestf
                bestf = f
                bestx .= x
            end

        end

        ind = i - pliminf + 1
        
        v[:, ind] .= bestx
        vf[ind]    = bestf

        println("Finished. p = $(i) and f = $(bestf).-> $(bestx)")
        
    end

    votsis = zeros(nInfo)
    
    for i = 1:nInfo
        for j = 1:nInfo
            if norm(v[:, i] - v[:, j]) < 10.0^(-3)
                votsis[i] += 1
            end
        end
    end
    
    mainind = findlast(x->x == maximum(votsis), votsis)

    println(v[:, mainind], ", ", vf[mainind], ",", pliminf + mainind - 1)
    
    return v[:, mainind], vf[mainind], pliminf + mainind - 1
    
end


end
