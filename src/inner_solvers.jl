export lmlovo, gnlslovo,lmlovoacelerated,lmlovomodified 

function lmlovoacelerated(model::Function, gmodel!::Function, θ::Vector{Float64},
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
    ResFun! = let

        data_    = data
        model_   = model
        gmodel!_ = gmodel!
        
        (θ::Vector{Float64}, ind, r::Vector{Float64},
         rJ::Array{Float64, 2}) -> begin
         
             nj += 1
         
             residual_fg!(model_, gmodel!_, data_, θ, ind, r, rJ)
         
         end
        
    end
    ResFun_acelerated! = let

        data_    = data
        model_   = model
        
        (θ::Vector{Float64}, ind, r::Vector{Float64}) -> begin
         
             nj += 1
         
             residual_fg_acelerated!(model_, data_, θ, ind, r)
         
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
        
        try

            # TODO: use a smart version of the QR factorization
            #F = qr!(G)
            F = cholesky!(G, Val(true))
            
            ldiv!(d, F, -grad_lovo)
    
            θ .= θ .+ d

            ResFun_acelerated!(θ, ind_lovo, val_res)

            BLAS.gemv!('T', 1.0, jac_res, val_res, 0.0, grad_lovo)
            
            ldiv!(d, F, -grad_lovo)
            
        catch

            with_logger(lm_logger.x) do
            
                @warn "Failed to solve the linear system. Will try new point."
                d .= - 1.0 .* grad_lovo 
                θ .= rand(n)
                
            end
            
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
function lmlovoacelerated(model::Function, θ::Vector{Float64}, data::Array{Float64,2},
                n::Int, p::Int; kwargs...)

    # Define closures for derivative and initializations

    # 'x' is considered as global parameter for this function
    model_cl(θ) = model(x, θ)
    
    grad_model!(g, x_, θ) = begin

        global x = x_
        
        ForwardDiff.gradient!(g, model_cl, θ)

    end

    return lmlovoacelerated(model, grad_model!, θ, data, n, p; kwargs...)
    
end

lmlovoacelerated(model::Function, gmodel!::Function, data::Array{Float64,2}, n::Int, p::Int; kwargs...) =
    lmlovo(model, gmodel!, zeros(Float64, n), data, n, p; kwargs...)

lmlovoacelerated(model::Function, data::Array{Float64,2}, n::Int, p::Int; kwargs...) =
    lmlovo(model, zeros(Float64, n), data, n, p; kwargs...)



function lmlovomodified(model::Function, gmodel!::Function, θ::Vector{Float64},
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
    ResFun! = let

        data_    = data
        model_   = model
        gmodel!_ = gmodel!
        
        (θ::Vector{Float64}, ind, r::Vector{Float64},
         rJ::Array{Float64, 2}) -> begin
         
             nj += 1
         
             residual_fg!(model_, gmodel!_, data_, θ, ind, r, rJ)
         
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
    θb = Vector{Float64}(undef, n)
    d = Vector{Float64}(undef, n)
    d1 = Vector{Float64}(undef, n)
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
        
        try

            # TODO: use a smart version of the QR factorization
            #F = qr!(G)
            F = cholesky!(G, Val(true))
            
            ldiv!(d1, F, grad_lovo)
            d1 .*= -1.0
    
            θb .= θ .+ d1 
            ResFun!(θb, ind_lovo, val_res, jac_res)

            BLAS.gemv!('T', 1.0, jac_res, val_res, 0.0, grad_lovo)
            
            ldiv!(d, F, grad_lovo)

            d .*= -1.0 

            d .= d1 .+ d
            
        catch

            with_logger(lm_logger.x) do
            
                @warn "Failed to solve the linear system. Will try new point."
                d .= - 1.0 .* grad_lovo 
                θ .= rand(n)
                
            end
            
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


function lmlovomodified(model::Function, θ::Vector{Float64}, data::Array{Float64,2},
                n::Int, p::Int; kwargs...)

    # Define closures for derivative and initializations

    # 'x' is considered as global parameter for this function
    model_cl(θ) = model(x, θ)
    
    grad_model!(g, x_, θ) = begin

        global x = x_
        
        ForwardDiff.gradient!(g, model_cl, θ)

    end

    return lmlovomodified(model, grad_model!, θ, data, n, p; kwargs...)
    
end

lmlovomodified(model::Function, gmodel!::Function, data::Array{Float64,2}, n::Int, p::Int; kwargs...) =
    lmlovo(model, gmodel!, zeros(Float64, n), data, n, p; kwargs...)

lmlovomodified(model::Function, data::Array{Float64,2}, n::Int, p::Int; kwargs...) =
    lmlovo(model, zeros(Float64, n), data, n, p; kwargs...)



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
    ResFun! = let

        data_    = data
        model_   = model
        gmodel!_ = gmodel!
        
        (θ::Vector{Float64}, ind, r::Vector{Float64},
         rJ::Array{Float64, 2}) -> begin
         
             nj += 1
         
             residual_fg!(model_, gmodel!_, data_, θ, ind, r, rJ)
         
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
        
        try

            # TODO: use a smart version of the QR factorization
            #F = qr!(G)
            F = cholesky!(G, Val(true))
            
            ldiv!(d, F, grad_lovo)

            d .*= -1.0
            
        catch

            with_logger(lm_logger.x) do
            
                @warn "Failed to solve the linear system. Will try new point."
                d .= - 1.0 .* grad_lovo 
                θ .= rand(n)
                
            end
            
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

    gnlslovo(model, gmodel!, θ, data::Array{T, 2}, n, p;
             ε::Number=1.0e-4, MAXITER=400, αls=2.0, dinc=2.0,
             MAXLSITER=100) where {T<:Float64}

    gnlslovo(model, θ::Vector{Float64}, data::Array{Float64,2},
             n::Int, p::Int; kwargs...)

    gnlslovo(model, gmodel!, data::Array{Float64,2}, n::Int,
             p::Int; kwargs...)

    gnlslovo(model, data::Array{Float64,2}, n::Int, p::Int; kwargs...)

LOVO Gauss-Newton with line-search described in

> R. Andreani, G. Cesar, R. M. Cesar-Jr., J. M. Martínez, and
> P. J. S. Silva, “Efficient curve detection using a {Gauss-Newton}
> method with applications in agriculture,” in Proc. 1st International
> Workshop on Computer Vision Applications for Developing Regions in
> Conjunction with ICCV 2007-CVDR-ICCV07, 2007.

Fit the `n`-parameter model `model` to the data given by matrix
`data`. The strategy is based on the LOVO function, which means that
only `p` (0 < `p` <= rows of `data`) points are trusted.

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
  - `αls`: number >1 to increase/decrease the parameter `t` in
    line-search
  - `dinc`: number >1 to increase the diagonal of the J^T J
    matrix in order to escape from singularity
  - `MAXLSITER`: maximum number of Linear System increases in
    diagonal before exiting. Also defines the maximum number of
    Line Search trials to satisfy Armijo (but does not exit in such
    case)

Returns a [`RAFFOutput`](@ref) object.

"""
function gnlslovo(model, gmodel!, θ, data::Array{T, 2}, n, p;
                  ε::Number=1.0e-4, MAXITER=400, αls=2.0, dinc=2.0,
                  MAXLSITER=100) where {T<:Float64}

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
    rp::Vector{Float64} = Vector{Float64}(undef, p)
    
    Jrp::Array{Float64, 2} = Array{Float64}(undef, p, n)
    
    # This function returns the residue and Jacobian of residue
    ResFun! = let

        data_    = data
        model_   = model
        gmodel!_ = gmodel!
        
        (θ::Vector{Float64}, ind, r::Vector{Float64},
         rJ::Array{Float64, 2}) -> begin
         
             nj += 1
         
             residual_fg!(model_, gmodel!_, data_, θ, ind, r, rJ)
         
         end
        
    end

    # Gauss-Newton algorithm

    maxoutind = min(p, 5)
    dtnf      = Inf
    status    = 1
    
    # Allocation
    θnew = Vector{Float64}(undef, n)
    d    = Vector{Float64}(undef, n)
    y    = Vector{Float64}(undef, n)
    G    = Array{Float64, 2}(undef, n, n)
    ∇f  = Vector{Float64}(undef, n)

    # Initial information
    
    ip, fk = LovoFun(θ)
    ResFun!(θ, ip, rp, Jrp)

    BLAS.gemv!('T', 1.0, Jrp, rp, 0.0, ∇f)

    norm∇f = norm(∇f, 2)

    safecount = 1

    safelscnt = 0

    # Main loop

    while (norm∇f >= ε) && (safecount < MAXITER)

        with_logger(lm_logger.x) do

            @info("Iteration $(safecount)")
            @info("  Current value:   $(fk)")
            @info("  ||∇f||_2: $(norm∇f)")
            @info("  Current iterate: $(θ)")
            @info("  Best indices (first $(maxoutind)): $(ip[1:maxoutind])")

        end

        # Solve the system

        BLAS.gemm!('T', 'N', 1.0, Jrp, Jrp, 0.0, G)

        safelscnt = 0
        
        while true && (safelscnt < MAXLSITER)
        
            F = qr(G)

            serror = false
            dtnf   = 0.0
            try
                ldiv!(d, F, ∇f)
                dtnf = - dot(d, ∇f)
            catch
                serror = true
            end

            with_logger(lm_logger.x) do
                @debug("  Found solution for linear system", d, dtnf)
            end

            (!serror) && (abs(dtnf) >= 1.0e-8) && break

            with_logger(lm_logger.x) do
                @debug("  Unable to solve the GN system (dtnf = $(dtnf))." *
                       "  Will increase the diagonal elements.")
            end

            # Adjust the diagonal elements of G in case of failure
            for i in diagind(G)
                G[i] = max(dinc * G[i], 1.0)
            end

            safelscnt += 1

        end

        if safelscnt >= MAXLSITER
            with_logger(lm_logger.x) do
                @debug("Unable to compute a descent direction.")
            end
            
            break
        end

        d .*= -1.0

        # Line search
        θnew .= θ .+ d
        
        t = 1.0

        safelscnt = 0

        # Warning. We cannot store ip to save computation, due to the
        # way that LovoFun is implemented (ip is a pointer to a hidden
        # vector, actually). This will change in the future
        ip, fkpd = LovoFun(θnew)

        while (fkpd <= fk + t * dtnf) && (safelscnt < MAXLSITER)
            t         *= αls
            θnew      .= θ .+ t .* d
            ip, fkpf   = LovoFun(θnew)
            safelscnt += 1
        end

        if t > 1.0
            t       /= αls
            θnew    .= θ .+ t .* d
            ip, fkpf = LovoFun(θnew)
        end

        # TODO: implement quadratic interpolation
        while (fkpd > fk + t * dtnf) && (safelscnt < MAXLSITER)
            t         /= αls
            θnew      .= θ .+ t .* d
            ip, fkpd   = LovoFun(θnew)
            safelscnt += 1
        end

        with_logger(lm_logger.x) do

            @info("""

                \tLine search
                \t  t          = $(t)
                \t  f(θ + t d) = $(fkpd)
                \t     d^T ∇f = $(dtnf)
                """)

            if safelscnt == MAXLSITER
                @warn("Armijo condition was not satisfied within" *
                      " $(MAXLSITER) iterations.")
            end

        end

        # New iteration

        θ .= θnew
        fk = fkpd

        ResFun!(θ, ip, rp, Jrp)

        BLAS.gemv!('T', 1.0, Jrp, rp, 0.0, ∇f)

        norm∇f = norm(∇f, 2)

        safecount += 1

    end

    if (norm∇f >= ε) && (safelscnt >= MAXLSITER)
        with_logger(lm_logger.x) do

            @info("Unable to compute a descent solution in " *
                  "($MAXLSITER) attempts.")

        end
        status = 0
    end        

    if (norm∇f >= ε) && (safecount == MAXITER)
        with_logger(lm_logger.x) do

            @info("No solution was found in $(safecount) iterations.")

        end
        status = 0
    end

    # TODO: Create a test for this case
    if isnan(norm∇f)
        with_logger(lm_logger.x) do

            @info("Incorrect value for gradient norm $(norm∇f).")

        end
        status = 0
    end

    safecount -= 1
    
    outliers = [1:npun;]
    setdiff!(outliers, ip)

    with_logger(lm_logger.x) do

        @info("""

        Final iteration (STATUS=$(status))
          Solution found:       $(θ)
          ||grad_lovo||_2:      $(norm∇f)
          Function value:       $(fk)
          Number of iterations: $(safecount)
          Outliers:             $(outliers)
    
        """)

    end
    
    return RAFFOutput(status, θ, safecount, p, fk, nf, nj, outliers)

end

function gnlslovo(model, θ::Vector{Float64}, data::Array{Float64,2},
                  n::Int, p::Int; kwargs...)

    # Define closures for derivative and initializations

    # 'x' is considered as global parameter for this function
    model_cl(θ) = model(x, θ)
    
    grad_model!(g, x_, θ) = begin

        global x = x_
        
        ForwardDiff.gradient!(g, model_cl, θ)

    end

    return gnlslovo(model, grad_model!, θ, data, n, p; kwargs...)
    
end

gnlslovo(model, gmodel!, data::Array{Float64,2}, n::Int, p::Int; kwargs...) =
    gnlslovo(model, gmodel!, zeros(Float64, n), data, n, p; kwargs...)

gnlslovo(model, data::Array{Float64,2}, n::Int, p::Int; kwargs...) =
    gnlslovo(model, zeros(Float64, n), data, n, p; kwargs...)

