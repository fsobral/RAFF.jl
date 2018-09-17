"""

This function is an auxiliary function. It finds the `p` smallest
values of vector `V` and brings them to the first `p` positions. The
indexes associated with the `p` smallest values are stored in `ind`.

"""
function SortFun!(V::Vector{Float64}, ind::Vector{Int}, p::Int)

    # If p is invalid, returns an empty view
    (p <= 0) && (return @view(ind[1:p]), @view(V[1:p]))
    
    npun = length(ind)
    
    for i = 1:npun
        ind[i] = i
    end
    
    for i = 1:p
        
        for j = i + 1:npun
            
            if V[i] > V[j]
                
                aux = ind[j]
                ind[j] = ind[i]
                ind[i] = aux
                
                vaux = V[j]
                V[j] = V[i]
                V[i] = vaux
            end
        end
    end

    return @view(ind[1:p]), @view(V[1:p])
    
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

"""

    generateNoisyData(model::Function, n::Int, np::Int, p::Int, xSol::Vector{Float64}=10.0 * rand(n),
                      tMin::Float64=-10.0, tMax::Float64=10.0)

Random generate a fitting one-dimensional data problem.

This function receives a `model(x, t)` function, the number of parameters
`n`, the number of points `np` to be generated and the number of
trusted points `p`. 

If the `n`-dimensional vector `xSol` is provided, the the exact
solution will not be random generated. The interval `[tMin, tMax]` for
generating the values to evaluate `model` can also be provided.

It returns a tuple `(data, xSol)` where

  - `data`: (`np` x `3`) array, where each row contains `t`,
    `model(xSol, t)` and `noise`.

  - `xSol`: `n`-dimensional vector with the exact solution.

"""
function generateNoisyData(model::Function, n, np, p, xSol::Vector{Float64}=10.0 * rand(n),
                           tMin::Float64=-10.0, tMax::Float64=10.0)

    @assert(tMin <= tMax, "Invalid interval for random number generation")
    
    #
    # Generate (ti,yi) where tMin <= t_i <= tMax (data)
    t = [tMin:(tMax - tMin) / (np - 1):tMax;]
    
    data = Array{Float64}(undef, np, 3)

    v = rand([1:1:np;], np - p)

    # Add noise to some random points
    for k = 1:np
            
        y = model(xSol, t[k])

        noise = 0.0
            
        if k in v 
            noise = randn() * rand([1.0, 2.0, 3.0])
        end
            
        data[k, 1] = t[k]
        data[k, 2] = y + noise
        data[k, 3] = noise

    end

    return data, xSol

end
