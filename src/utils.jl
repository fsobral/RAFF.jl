export generateTestProblems

"""

    eliminate_local_min!(sols::Vector{RAFFOutput})

Check if the function value of the solution found by smaller values of
`p` is not greater when compared with larger ones. This certainly
indicates that a local minimizer was found by the smaller `p`.

"""
function eliminate_local_min!(model::Function, data::Array{Float64, 2},
                              sols::Vector{RAFFOutput})

    lv = length(sols)

    sec_sol_ind = -1

    # Start from the largest p
    for i = lv:-1:1

        (sols[i].status != 1) && continue
        
        for j = i - 1:-1:1

            if (sols[j].status == 1) && (sols[j].f > sols[i].f)

                @debug("  Possible local minimizer for p = $(sols[j].p) " *
                       "with f = $(sols[j].f). Removing it.")
                
                sols[j] = RAFFOutput(0)

            end

        end

    end

    # Test the maximum 'p'

    nump, = size(data)
    
    maxp = sols[lv]
        
    if (lv > 1) && (maxp.status == 1)

        i = lv - 1

        # while (i > 0) && (sols[i].status != 1)

        #     i -= 1

        # end

        bestf = maxp.f
        j = 0
        
        while (i > 0)

            if  (sols[i].status == 1) && (sols[i].f < bestf)

                bestf = sols[i].f

                j = i

            end

            i -= 1

        end

        if j == 0
            
            @debug("  Nobody to reject p = $(maxp.p).")

        else

            sec_sol = sols[j]
            
            @debug("  p = $(sec_sol.p) will try to reject p = $(maxp.p).")

            nmin = 0

            for i = 1:nump

                y  = data[i, 2]
                y1 = model(maxp.solution, data[i, 1])
                y2 = model(sec_sol.solution, data[i, 1])
                
                (abs(y - y1) < abs(y - y2)) && (nmin += 1)

            end

            if nmin < lv / 2

                @debug("  nmin = $(nmin). Rejecting p = $(maxp.p).")

                sols[lv] = RAFFOutput(0)

            end

        end
        
    end
    
end


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
                              tMin=-10.0, tMax=10.0,
                              xSol=10.0 * randn(n), std=200.0, outTimes=7.0)
    
    open(solFilename, "w") do sol
    
        println(sol, n) # number of variables

        println(sol, xSol) # parameters
        
        println(sol, modelStr) # function expression

    end

    # Choose exactly np - p random points to be perturbed
    v = Vector{Int}(undef, np - p)

    ntp = 0

    points = [1:np;]

    while ntp < np - p
        
        u = rand(points, np - p - ntp)

        unique!(u)

        for i in u
            v[ntp + 1] = i
            ntp += 1
        end

        setdiff!(points, u)
        
    end

    #
    # Generate (ti,yi) where tMin <= t_i <= tMax (data)
    t = [tMin:(tMax - tMin) / (np - 1):tMax;]
    
    data = open(datFilename, "w")

    # Dimension of the domain of the function to fit
    @printf(data, "%d\n", 1)
    
    # Add noise to some random points
    for k = 1:np
            
        y = model(xSol, t[k]) + randn() * std

        noise = 0.0
            
        if k in v 
            y = model(xSol, t[k])
            noise = outTimes * std * sign(randn())
        end
            
        @printf(data, "%20.15f %20.15f %1d\n",
                t[k], y + noise, Int(noise != 0))

    end

    close(data)
    
end

"""

    generateNoisyData(model::Function, n::Int, np::Int, p::Int, xSol::Vector{Float64}=10.0 * rand(n),
                      tMin::Float64=-10.0, tMax::Float64=10.0)

    generateNoisyData(model::Function, n, np, p, tMin::Float64, tMax::Float64)

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
function generateNoisyData(model::Function, n::Int, np::Int, p::Int,
                           xSol::Vector{Float64}=10.0 * rand(n),
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

generateNoisyData(model::Function, n::Int, np::Int, p::Int, tMin::Float64, tMax::Float64) =
    generateNoisyData(model, n, np, p, 10.0 * rand(n), tMin, tMax)
