export lmlovo

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
                MAXITER::Int=400, ε::Float64=1.0e-4)

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

