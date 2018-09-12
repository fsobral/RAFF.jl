module RAFF

__precompile__(false)

# Dependencies
using ForwardDiff
using LinearAlgebra
using Printf

export LMlovo, raff, generateTestProblems

"""
    LMlovo(model::Function, data::Array{Float64, 2}, n::Int, p::Int [; kwargs...])

    LMlovo(model::Function, gmodel!::Function, data::Array{Float64,2}, n::Int,
           p::Int [; MAXITER::Int])

Fit the `n`-parameter model `model` to the data given by matrix
`data`. The strategy is based on the LOVO function, which means that
only `p` (0 < `p` <= `n`) points are trusted. The Levenberg-Marquardt
algorithm is implemented in this version.

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

"""
function LMlovo(model::Function, data::Array{Float64,2}, n::Int, p::Int;
                MAXITER::Int=200)

    # Define closures for derivative and initializations

    # 't' is considered as global parameter for this function
    model_cl(x) = model(x, t)
    
    grad_model(x, t_, g) = begin

        global t = t_
        
        ForwardDiff.gradient!(g, model_cl, x)

    end

    return LMlovo(model, grad_model, data, n, p, MAXITER=MAXITER)
    
end

function LMlovo(model::Function, gmodel!::Function,
                data::Array{Float64,2}, n::Int, p::Int;
                MAXITER::Int=200)

    @assert(n > 0, "Dimension should be positive.")
    @assert(p >= 0, "Trusted points should be nonnegative.")
    
    t = 0.0
    npun = length(data[:, 1])
    
    #All evaluated instance depends of a sort function 
    SortFun(V::Vector) = begin
        aux=0
        vaux=0.0
        ind=[1:1:npun;]
        for i=1:p
            for j=i+1:npun
                if (V[i]>V[j])
                    aux=ind[j]
                    ind[j]=ind[i]
                    ind[i]=aux
                    vaux=V[j]
                    V[j]=V[i]
                    V[i]=vaux
                end
            end
        end
        return ind[1:1:p]
    end
    
    # Main function - the LOVO function
    LovoFun(x) = begin #return a ordered set index and lovo value 
        F = zeros(npun)
        for i = 1:npun 
            F[i] = (model(x, data[i,1]) - data[i, 2])^2
        end
        return SortFun(F), sum(F[SortFun(F)])
    end
    
    # This function returns the residue and Jacobian of residue 
    ResFun(x, ind) = begin
        r = zeros(p)
        rJ = zeros(p, n)
        k = 1
        for i in ind
            t = data[i, 1]
            r[k] = model(x, t) - data[i, 2]

            v = @view(rJ[k, :])
            gmodel!(x, t, v)

            k = k + 1
        end
        return r, rJ
    end
    
    # Levenberg-Marquardt algorithm
    # -----------------------------
    
    Id = Matrix(1.0I, n, n)
    
    # Parameters
    ε = 10.0^(-4)
    λ_up = 2.0
    λ_down = 2.0
    λ = 1.0

    x = zeros(n) #initial point
    d = zeros(n)
    y = zeros(n)
    (ind_lovo,val_lovo)=LovoFun(x)
    (val_res,jac_res)=ResFun(x,ind_lovo)
    grad_lovo=jac_res'*val_res
    safecount=1
    while (norm(grad_lovo,2) >= ε) && (safecount < MAXITER)
        G = jac_res' * jac_res + λ * Id
        F = qr(G)     
        ad = try
            y = F.Q \ (- grad_lovo)
            d = F.R \ y
        catch
            "error"
        end
        if ad == "error" #restarting if lapack fails
            d = - grad_lovo 
            x = rand(n)
        else 
            d = ad
        end
        #d=G\(-grad_lovo)
        xnew=x+d
        (ind_lovo_new,val_lovo_new)=LovoFun(xnew)
        if  val_lovo_new<=val_lovo
            x=copy(xnew)
            val_lovo=copy(val_lovo_new)
            ind_lovo=copy(ind_lovo_new)
            λ=λ/(λ_down)
            (val_res,jac_res)=ResFun(x,ind_lovo)
            grad_lovo=jac_res'*val_res
        else
            λ=λ*(λ_up)
        end
        safecount+=1
    end
    if safecount == MAXITER
    #    println("no solution was found in $safecount iterations")
        return 0, x, safecount, p
    else
    #    println("solution found::   $x " )
    #    println("number of iterations:: $(safecount)")
        return 1, x, safecount, p
    end

end

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
    
    k = 1

    for i = pliminf:plimsup
        v[k] = LMlovo(model, gmodel!, data, n, i)
        k += 1
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

    generateTestProblems(datFilename::String, solFilename::String,
                         model::Function, modelStr::String, n::Int,
                         np::Int, p::Int)

Generate random data files for testing fitting problems.

  - `datFilename` and `solFilename` are strings with the name of the
    files for storing the random data and solution, respectively.
  - `model` is the model function and `modelStr` is a string
    representing this model function, e.g.

         model = (x, t) -> x[1] * t + x[2]
         modelStr = "(x, t) -> x[1] * t + x[2]"

    where `x` represents the parameters (to be found) of the model and
    `t` is the variable of the model.
  - `n` is the number of parameters
  - `np` is the number of points to be generated.
  - `p` is the number of trusted points to be used in the LOVO
    approach.

"""
function generateTestProblems(datFilename::String,
                              solFilename::String, model::Function,
                              modelStr::String,
                              n::Int, np::Int, p::Int;
                              tmin=-10.0, tmax=10.0)
                           
    # Generate parameters x (solution)
    x = 10.0 * randn(n)
    
    open(solFilename, "w") do sol
    
        println(sol, n) # number of variables

        println(sol, x) # parameters
        
        println(sol, modelStr) # function expression

    end

    #
    # Generate (ti,yi) where tmin <= t_i <= tmax (data)
    t = [tmin:(tmax - tmin) / (np - 1):tmax;]
    
    data = open(datFilename, "w")

    v = rand([1:1:np;], np - p)

    # Add noise to some random points
    for k = 1:np
            
        y = model(x, t[k])

        noise = 0.0
            
        if k in v 
            noise = randn() * rand([1.0, 2.0, 3.0])
        end
            
        @printf(data, "%20.15f %20.15f %20.15f\n",
                t[k], y + noise, noise)

    end

    close(data)
    
end


end
