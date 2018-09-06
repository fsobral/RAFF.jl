module RAFF

# Dependencies
using ForwardDiff
using LinearAlgebra

export LMlovo, raff

"""
    LMlovo(model, data, n, p)

Fit the `n`-parameter model `model` to the data given by matrix
`data`. The strategy is based on the LOVO function, which means that
only `p` (0 < `p` <= `n`) points are trusted. The Levenberg-Marquardt
algorithm is implemented in this version.

The signature of function `model` should be given by

    model(x, t)

where `x` is a `n`-dimensional vector of parameters and `t` is the
argument.

Returns a tuple `s`, `x`, `iter`, `p`, where

  - `s`: is 1 if converged and 0 if not
  - `x`: vector with the parameters of the model
  - `iter`: number of iterations up to convergence
  - `p`: number of trusted points

"""
function LMlovo(model::Function, data::Array{Float64,2}, n::Int, p::Int)

    # Define closures for derivative and initializations
    model_cl(x) = model(x, t)
    grad_model_cl(x) = ForwardDiff.gradient(model_cl, x)
    
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
            r[k] = model(x, data[i, 1]) - data[i, 2]
            t = data[i, 1]
            rJ[k, :] = grad_model_cl(x)
            k = k + 1
        end
        return r, rJ
    end
    # Levenberg-Marquardt algorithm
    Id = Matrix(1.0I, n, n)
    ε=10.0^(-4)
    λ_up=2.0
    λ_down=2.0
    λ=1.0
    x=zeros(n) #initial point
    d=zeros(n)
    y=zeros(n)
    (ind_lovo,val_lovo)=LovoFun(x)
    (val_res,jac_res)=ResFun(x,ind_lovo)
    grad_lovo=jac_res'*val_res
    safecount=1
    while norm(grad_lovo,2)>=ε && safecount<200
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
    if safecount==200
    #    println("no solution was founded in $safecount iterations")
        return 0,x,safecount,p
    else
    #    println("solution founded::   $x " )
    #    println("number of iterations:: $(safecount)")
        return 1,x,safecount,p
    end
end

"""
    raff(model, data, n)

  - `model`: function to fit data. Its signature should be given by

      model(x, t)

    where `x` is a `n`-dimensional vector of parameters and `t` is the
    argument
  - `data`: data to be fit
  - `n`: dimension of the parameter vector in the model function

Robust Algebric Fitting Function (RAFF) algorithm. This function uses
a voting system to automatically find the number of trusted data
points to fit the `model`.

"""
function raff(model, data, n)

    pliminf = Int(round(length(data[:, 1]) / 2.0))
    plimsup = length(data[:, 1])
    v = Array{Any,1}(undef, plimsup - pliminf + 1)
    k = 1

    for i = pliminf:plimsup
        v[k] = LMlovo(model, data, n, i)
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
    
    display(v[mainind])
end

end
