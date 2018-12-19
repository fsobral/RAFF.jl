export generateTestProblems, generateNoisyData

# This dictionary represents the list of models used in the tests
# Return the tuple (n, model, model_str)
const model_list = Dict(
    "linear" => (2, (x, t) -> x[1] * t[1] + x[2],
                 "(x, t) -> x[1] * t[1] + x[2]"),
    "cubic" => (4, (x, t) -> x[1] * t[1]^3 + x[2] * t[1]^2 + x[3] * t[1] + x[4],
                "(x, t) -> x[1] * t[1]^3 + x[2] * t[1]^2 + x[3] * t[1] + x[4]"),
    "expon" => (3, (x, t) -> x[1] + x[2] * exp(- x[3] * t[1]),
                "(x, t) -> x[1] + x[2] * exp(- x[3] * t[1])"),
    "logistic" => (4, (x, t) -> x[1] + x[2] / (1.0 + exp(- x[3] * t[1] + x[4])),
                   "(x, t) -> x[1] + x[2] / (1.0 + exp(- x[3] * t[1] + x[4]))")
)


"""

    generateTestProblems(datFilename::String, solFilename::String,
                         model::Function, modelStr::String, n::Int,
                         np::Int, p::Int)

Generate random data files for testing fitting problems.

  - `datFilename` and `solFilename` are strings with the name of the
    files for storing the random data and solution, respectively.
  - `model` is the model function and `modelStr` is a string
    representing this model function, e.g.

         model = (x, t) -> x[1] * t[1] + x[2]
         modelStr = "(x, t) -> x[1] * t[1] + x[2]"

    where vector `x` represents the parameters (to be found) of the
    model and vector `t` are the variables of the model.
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

    # Generate solution file
    
    open(solFilename, "w") do sol
    
        println(sol, n) # number of variables

        println(sol, xSol) # parameters
        
        println(sol, modelStr) # function expression

    end

    # Generate data file
    
    open(datFilename, "w") do data
    
        vdata, xsol, outliers = generateNoisyData(model, n, np, p;
                                   tMin=tMin, tMax=tMax, xSol=xSol,
                                   std=std, outTimes=outTimes)
    
        # Dimension of the domain of the function to fit
        @printf(data, "%d\n", 1)

        for k = 1:np

            @printf(data, "%20.15f %20.15f %1d\n",
                    vdata[k, 1], vdata[k, 2], Int(k in outliers))

        end

    end
    
end

"""

    get_unique_random_points(np::Int, npp::Int)

Choose exactly `npp` unique random points from a set containing `np`
points. This function is similar to `rand(vector)`, but does not allow
repetitions.

Return a vector with the selected points.

"""
function get_unique_random_points(np::Int, npp::Int)

    # Check invalid arguments
    ((np <= 0) || (npp <=0)) && (return Vector{Int}())

    ntp = min(npp, np)
    
    v = Vector{Int}(undef, ntp)

    points = [1:np;]

    while ntp > 0
        
        u = rand(points, ntp)

        unique!(u)

        for i in u
            v[ntp] = i
            ntp -= 1
        end

        setdiff!(points, u)
        
    end

    return v

end

"""

    generateNoisyData(model::Function, n::Int, np::Int, p::Int;
                      tMin::Float64=-10.0, tMax::Float64=10.0,
                      xSol::Vector{Float64}=10.0 * randn(Float64, n),
                      std::Float64=200.0, outTimes::Float64=7.0)

    generateNoisyData(model::Function, n, np, p, tMin::Float64, tMax::Float64)

    generateNoisyData(model::Function, n::Int, np::Int, p::Int,
                      xSol::Vector{Float64}, tMin::Float64, tMax::Float64)

Random generate a fitting one-dimensional data problem.

This function receives a `model(x, t)` function, the number of parameters
`n`, the number of points `np` to be generated and the number of
trusted points `p`. 

If the `n`-dimensional vector `xSol` is provided, the the exact
solution will not be random generated. The interval `[tMin, tMax]` for
generating the values to evaluate `model` can also be provided.

It returns a tuple `(data, xSol, outliers)` where

  - `data`: (`np` x `2`) array, where each row contains `t` and
    `model(xSol, t)`.
  - `xSol`: `n`-dimensional vector with the exact solution.
  - `outliers`: the outliers of this data set

"""
function generateNoisyData(model::Function, n::Int, np::Int, p::Int;
                           tMin::Float64=-10.0, tMax::Float64=10.0,
                           xSol::Vector{Float64}=10.0 * randn(Float64, n),
                           std::Float64=200.0, outTimes::Float64=7.0)

    @assert(tMin <= tMax, "Invalid interval for random number generation")
    
    # Generate (ti,yi) where tMin <= t_i <= tMax (data)
    t = [tMin:(tMax - tMin) / (np - 1):tMax;]
    
    data = Array{Float64}(undef, np, 3)

    # Points selected to be outliers
    v = get_unique_random_points(np, np - p)
    
    # Add noise to some random points
    for k = 1:np
            
        y = model(xSol, t[k]) + randn() * std

        noise = 0.0
            
        if k in v 
            y = model(xSol, t[k])
            noise = outTimes * std * sign(randn())
        end
            
        data[k, 1] = t[k]
        data[k, 2] = y + noise
        data[k, 3] = noise

    end

    return data, xSol, v

end

generateNoisyData(model::Function, n::Int, np::Int, p::Int, tMin::Float64, tMax::Float64) =
    generateNoisyData(model, n, np, p; tMin=tMin, tMax=tMax)

generateNoisyData(model::Function, n::Int, np::Int, p::Int,
                  xSol::Vector{Float64}, tMin::Float64, tMax::Float64) =
                      generateNoisyData(model, n, np, p; tMin=tMin, tMax=tMax, xSol=xSol)
