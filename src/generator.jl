export generate_test_problems, generate_noisy_data!,
    generate_noisy_data, generate_clustered_noisy_data,
    generate_clustered_noisy_data!, generate_circle,
    generate_ncircle

"""

This dictionary represents the list of models used in the generation of random tests.
Return the tuple `(n, model, model_str)`, where

  - `n` is the number of parameters of the model
  - `model` is the model of the form `m(x, θ)`, where `x` are the
    variables and `θ` are the parameters
  - `model_str` is the string representing the model, used to build random generated problems

"""
const model_list = Dict(
    "linear" => (2, (x, θ) -> θ[1] * x[1] + θ[2],
                 "(x, θ) -> θ[1] * x[1] + θ[2]"),
    "cubic" => (4, (x, θ) -> θ[1] * x[1]^3 + θ[2] * x[1]^2 + θ[3] * x[1] + θ[4],
                "(x, θ) -> θ[1] * x[1]^3 + θ[2] * x[1]^2 + θ[3] * x[1] + θ[4]"),
    "expon" => (3, (x, θ) -> θ[1] + θ[2] * exp(- θ[3] * x[1]),
                "(x, θ) -> θ[1] + θ[2] * exp(- θ[3] * x[1])"),
    "logistic" => (4, (x, θ) -> θ[1] + θ[2] / (1.0 + exp(- θ[3] * x[1] + θ[4])),
                   "(x, θ) -> θ[1] + θ[2] / (1.0 + exp(- θ[3] * x[1] + θ[4]))"),
    "circle" => (3, (x, θ) -> (x[1] - θ[1])^2 + (x[2] - θ[2])^2 - θ[3]^2,
                 "(x, θ) -> (x[1] - θ[1])^2 + (x[2] - θ[2])^2 - θ[3]^2"),
    "ellipse" => (6, (x, θ) -> θ[1] * x[1]^2 + θ[2] * x[1] * x[2] + θ[3] * x[2]^2 +
                  θ[4] * x[1] + θ[5] * x[2] + θ[6],
                  "(x, θ) -> θ[1] * x[1]^2 + θ[2] * x[1] * x[2] + θ[3] * x[2]^2 + θ[4] * x[1] + θ[5] * x[2] + θ[6]")
)

"""

    interval_rand!(x::Vector{Float64},
        intervals::Vector{Tuple{Float64, Float64}})

Fill a vector `x` with uniformly distributed random numbers generated
in the interval given by `intervals`. It is assumed that `length(x) ==
length(intervals)`.

Throws an `ErrorException` if the dimension of `x` is smaller the
dimension of `intervals` or if the intervals are invalid.

"""
function interval_rand!(x::Vector{Float64},
    intervals::Vector{Tuple{Float64, Float64}})

    (length(x) < length(intervals)) &&
        error("Length of vector smaller than length of intervals.")

    for v in intervals

        (v[1] > v[2]) && error("Bad interval $(v[1])>$(v[2]).")

    end
    
    map!((v) -> v[1] + rand() * (v[2] - v[1]), x, intervals)
    
end

"""

    generate_test_problems(datFilename::String, solFilename::String,
        model::Function, modelStr::String, n::Int, np::Int, p::Int;
        x_interval::Tuple{Float64, Float64}=(-10.0, 10.0),
        θSol::Vector{Float64}=10.0 * randn(n), std::Float64=200.0,
        out_times::Float64=7.0)

    generate_test_problems(datFilename::String, solFilename::String,
        model::Function, modelStr::String, n::Int, np::Int, p::Int,
        cluster_interval::Tuple{Float64, Float64};
        x_interval::Tuple{Float64, Float64}=(-10.0, 10.0),
        θSol::Vector{Float64}=10.0 * randn(n), std::Float64=200.0,
        out_times::Float64=7.0)

Generate random data files for testing fitting problems.

  - `datFilename` and `solFilename` are strings with the name of the
    files for storing the random data and solution, respectively.
  - `model` is the model function and `modelStr` is a string
    representing this model function, e.g.

         model = (x, θ) -> θ[1] * x[1] + θ[2]
         modelStr = "(x, θ) -> θ[1] * x[1] + θ[2]"

    where vector `θ` represents the parameters (to be found) of the
    model and vector `x` are the variables of the model.
  - `n` is the number of parameters
  - `np` is the number of points to be generated.
  - `p` is the number of trusted points to be used in the LOVO
    approach.

If `cluster_interval` is provided, then generates outliers only in
this interval.

Additional parameters:

  - `xMin`, `xMax`: interval for generating points in one dimensional
    tests *Deprecated*
  - `x_interval`: interval for generating points in one dimensional
    tests
  - `θSol`: true solution, used for generating perturbed points
  - `std`: standard deviation
  - `out_times`: deviation for outliers will be `out_times * std`.

"""
function generate_test_problems(datFilename::String,
    solFilename::String, model::Function, modelStr::String, n::Int,
    np::Int, p::Int; gn_kwargs...)

    # Generate data file
    
    θSol = nothing
    
    open(datFilename, "w") do data
    
        vdata, θSol, outliers = generate_noisy_data(model, n, np, p;
                                                    gn_kwargs...)
    
        # Dimension of the domain of the function to fit
        @printf(data, "%d\n", 1)

        for k = 1:np

            @printf(data, "%20.15f %20.15f %1d\n",
                    vdata[k, 1], vdata[k, 2], Int(k in outliers))

        end

    end

    # Generate solution file
    
    open(solFilename, "w") do sol
    
        println(sol, n) # number of variables

        println(sol, θSol) # parameters
        
        println(sol, modelStr) # function expression

    end

end

function generate_test_problems(datFilename::String,
    solFilename::String, model::Function, modelStr::String, n::Int,
    np::Int, p::Int, x_interval::Tuple{Float64, Float64},
    cluster_interval::Tuple{Float64, Float64}; gn_kwargs...)

    # Generate data file
    
    θSol = nothing
    
    open(datFilename, "w") do data
    
        vdata, θSol, outliers = generate_clustered_noisy_data(model,
            n, np, p, x_interval, cluster_interval; gn_kwargs...)
    
        # Dimension of the domain of the function to fit
        @printf(data, "%d\n", 1)

        for k = 1:np

            @printf(data, "%20.15f %20.15f %1d\n",
                    vdata[k, 1], vdata[k, 2], Int(k in outliers))

        end

    end

    # Generate solution file
    
    open(solFilename, "w") do sol
    
        println(sol, n) # number of variables

        println(sol, θSol) # parameters
        
        println(sol, modelStr) # function expression

    end

end

"""

    generate_circle(dat_filename::String, np::Int, p::Int;
        std::Float64=0.1, θSol::Vector{Float64}=1.0*randn(Float64, 3),
        outTimes::Float64=3.0, interval=(rand(i)*2.0*π for i = 1:np))

Generate perturbed points in a circle given by `θSol` and save to
`dat_filename` in RAFF format. Return the np x 4 matrix with data (the
4th column is 0 if the point is "correct") and a `np - p` integer
vector containing the points selected to be outliers.

  - `dat_filename` is a String with the name of the file to store
    generated data.
  - `np` is the number of points to be generated.
  - `p` is the number of *trusted points* to be used in the LOVO
    approach.

Additional configuration parameters are

  - `std`: standard deviation.
  - `θSol`: true solution, used for generating perturbed points.
  - `out_times`: deviation for outliers will be `out_times * std`.
  - `interval`: any iterable object containing `np` numbers between 0
    and 2π.

"""
function generate_circle(dat_filename::String, np::Int, p::Int;
    std::Float64=0.1, θSol::Vector{Float64}=1.0*randn(Float64, 3),
    outTimes::Float64=3.0, interval=(rand()*2.0*π for i = 1:np))

    (length(interval) != np) &&
        error("Size of interval different from given value of np")

    ρ = (α, ρ) -> [ρ * cos(α) + θSol[1], ρ * sin(α) + θSol[2]]

    data = Array{Float64, 2}(undef, np, 4)

    # Points selected to be outliers
    v = RAFF.get_unique_random_points(np, np - p)

    for (i, α) in enumerate(interval)

        pt = ρ(α, θSol[3] + std * randn())

        data[i, 3:4] .= 0.0

        # Follow the noise idea of J. Yu, H. Zheng, S. R. Kulkarni,
        # and H. V. Poor, "Two-Stage Outlier Elimination for Robust
        # Curve and Surface Fitting," EURASIP J. Adv. Signal Process.,
        # vol. 2010, no. 1, p. 154891, Dec. 2010.
        if i in v

            pt[1] += outTimes * std * randn()
            pt[2] += outTimes * std * randn()

            data[i, 4] = 1.0

        end

        data[i, 1:2] = pt

    end

    open(dat_filename, "w") do fp
    
        # Dimension of the domain of the function to fit
        @printf(fp, "%d\n", 2)

        for k = 1:np

            @printf(fp, "%20.15f %20.15f %20.15f %1d\n",
                    data[k, 1], data[k, 2], data[k, 3], Int(k in v))

        end

    end

    return data, v

end

"""

    generate_ncircle(dat_filename::String,np::Int, p::Int;
      std::Float64=0.1, θSol::Vector{Float64}=10.0*randn(Float64, 3),
      interval=(rand()*2.0*π for i = 1:np))

Generate perturbed points and uniform noise in a square containing the
circle given by `θSol` and save data to `dat_filename` in RAFF
format. Return the np x 4 matrix with data (the 4th column is 0 if the
point is "correct") and a `np - p` integer vector containing the
points selected to be outliers.

  - `dat_filename` is a String with the name of the file to store
    generated data.
  - `np` is the number of points to be generated.
  - `p` is the number of *trusted points* to be used in the LOVO
    approach.

Additional configuration parameters are

  - `std`: standard deviation.
  - `θSol`: true solution, used for generating perturbed points.
  - `interval`: any iterable object containing `np` numbers between 0
    and 2π.

"""
function generate_ncircle(dat_filename::String, np::Int, p::Int;
    std::Float64=0.1, θSol::Vector{Float64}=10.0*randn(Float64, 3),
    interval=(rand()*2.0*π for i = 1:np))

    (length(interval) != np) &&
        error("Size of interval different from given value of np")

    ρ = (α, ρ) -> [ρ * cos(α) + θSol[1], ρ * sin(α) + θSol[2]]

    data = Array{Float64, 2}(undef, np, 4)

    for (i, α) in enumerate(interval)

        pt = ρ(α, θSol[3] + std * randn())

        data[i, 1:2]  = pt
        data[i, 3:4] .= 0.0

    end

    # Just random noise

    v = Vector{Int}(undef, np - p)

    for i = p + 1:np

        data[i, 1] = θSol[1] - 2.0 * θSol[3] + rand() * 4.0 * θSol[3]
        data[i, 2] = θSol[2] - 2.0 * θSol[3] + rand() * 4.0 * θSol[3]
        data[i, 3] = 0.0
        data[i, 4] = 1.0

        v[i - p] = i

    end

    open(dat_filename, "w") do fp

        # Dimension of the domain of the function to fit
        @printf(fp, "%d\n", 2)

        for k = 1:np

            @printf(fp, "%20.15f %20.15f %20.15f %1d\n",
                    data[k, 1], data[k, 2], data[k, 3], Int(k in v))

        end

    end

    return data, v

end

"""

    get_unique_random_points(np::Int, npp::Int)

Choose exactly `npp` unique random points from a set containing `np`
points. This function is similar to `rand(vector)`, but does not allow
repetitions.

If `npp` < `np`, returns all the `np` points. Note that this function
is not very memory efficient, since the process of selecting unique
elements involves creating several temporary vectors.

Return a vector with the selected points.

"""
function get_unique_random_points(np::Int, npp::Int)

    # Check invalid arguments
    ((np <= 0) || (npp <=0)) && (return Vector{Int}())

    ntp = min(npp, np)
    
    v = Vector{Int}(undef, ntp)
    
    return get_unique_random_points!(v, np, npp)

end

"""

    get_unique_random_points!(v::Vector{Int}, np::Int, npp::Int)

Choose exactly `npp` unique random points from a set containing `np`
points. This function is similar to `rand(vector)`, but does not allow
repetitions.

If `npp` < `np`, returns all the `np` points. Note that this function
is not very memory efficient, since the process of selecting unique
elements involves creating several temporary vectors.

Return the vector `v` provided as argument filled with the selected
points.

"""
function get_unique_random_points!(v::Vector{Int}, np::Int, npp::Int)

    # Check invalid arguments
    ((np <= 0) || (npp <=0)) && return v

    ntp = min(np, npp)
    
    # Check invalid size
    (length(v) < ntp) && throw(DimensionMismatch("Incorrect size for vector."))
    
    # Check small np for efficiency
    if np == ntp

        for i = 1:np
            v[i] = i
        end

        return v

    end

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

    generate_noisy_data(model::Function, n::Int, np::Int, p::Int;
        x_interval::Tuple{Float64, Float64}=(-10.0, 10.0),
        θSol::Vector{Float64}=10.0 * randn(Float64, n),
        std::Float64=200.0, out_times::Float64=7.0)

    generate_noisy_data(model::Function, n::Int, np::Int, p::Int,
        x_interval::Tuple{Float64, Float64})

    generate_noisy_data(model::Function, n::Int, np::Int, p::Int,
        θSol::Vector{Float64}, x_interval::Tuple{Float64, Float64})


Random generate a fitting one-dimensional data problem.

This function receives a `model(x, θ)` function, the number of parameters
`n`, the number of points `np` to be generated and the number of
trusted points `p`. 

If the `n`-dimensional vector `θSol` is provided, then the exact
solution will not be random generated. The interval `[xMin, xMax]`
(*deprecated*) or `x_interval` for generating the values to evaluate
`model` can also be provided.

It returns a tuple `(data, θSol, outliers)` where

  - `data`: (`np` x `3`) array, where each row contains `x` and
    `model(x, θSol)`.
  - `θSol`: `n`-dimensional vector with the exact solution.
  - `outliers`: the outliers of this data set

"""
function generate_noisy_data(model::Function, n::Int, np::Int,
    p::Int; gn_kwargs...)

    # Data matrix
    data = Array{Float64}(undef, np, 3)

    # Points selected to be outliers
    v = Vector{Int}(undef, np - p)

    return generate_noisy_data!(data, v, model, n, np, p;
                                gn_kwargs...)
    
end

# Deprecated

@deprecate(generate_noisy_data(model::Function, n, np, p,
    xMin::Float64, xMax::Float64), generate_noisy_data(model, n, np,
    p, (xMin, xMax)))

@deprecate(generate_noisy_data(model::Function, n::Int, np::Int,
    p::Int, θSol::Vector{Float64}, xMin::Float64, xMax::Float64),
    generate_noisy_data(model, n, np, p, θSol, (xMin, xMax)))

generate_noisy_data(model::Function, n::Int, np::Int, p::Int,
    x_interval::Tuple{Float64, Float64}) = generate_noisy_data(model,
    n, np, p; x_interval=x_interval)

generate_noisy_data(model::Function, n::Int, np::Int, p::Int,
    θSol::Vector{Float64}, x_interval::Tuple{Float64, Float64}) =
    generate_noisy_data(model, n, np, p; x_interval=x_interval,
    θSol=θSol)

"""

    generate_noisy_data!(data::AbstractArray{Float64, 2},
        v::Vector{Int}, model::Function, n::Int, np::Int, p::Int;
        x_interval::Tuple{Float64, Float64}=(-10.0, 10.0),
        θSol::Vector{Float64}=10.0 * randn(Float64, n),
        std::Float64=200.0, out_times::Float64=7.0)

Random generate a fitting one-dimensional data problem, storing the
data in matrix `data` and the outliers in vector `v`.

This function receives a `model(x, θ)` function, the number of parameters
`n`, the number of points `np` to be generated and the number of
trusted points `p`. 

If the `n`-dimensional vector `θSol` is provided, then the exact
solution will not be random generated. The interval `[xMin, xMax]`
(*deprecated*) or `x_interval` for generating the values to evaluate
`model` can also be provided.

It returns a tuple `(data, θSol, outliers)` where

  - `data`: (`np` x `3`) array, where each row contains `x` and
    `model(x, θSol)`.
  - `θSol`: `n`-dimensional vector with the exact solution.
  - `outliers`: the outliers of this data set

"""
function generate_noisy_data!(data::AbstractArray{Float64, 2},
    v::Vector{Int}, model::Function, n::Int, np::Int, p::Int;
    x_interval::Tuple{Float64, Float64}=(-10.0, 10.0),
    θSol::Vector{Float64}=10.0 * randn(Float64, n),
    std::Float64=200.0, out_times::Float64=7.0)

    @assert(x_interval[1] <= x_interval[2],
            "Invalid interval for random number generation.")
    @assert(size(data) == (np, 3),
            "Invalid size of data matrix. $(size(data)) != $((np, 3)).")
    @assert(length(v) >= np - p,
            "Invalid size for vector of outliers.")
    
    # Generate (x_i) where x_interval[1] <= x_i <= x_interval[2] (data)
    # Fix the problem of large interval with 1 element.
    x = (np == 1) ? sum(x_interval) / 2.0 : LinRange(x_interval[1], x_interval[2], np)
    
    # Points selected to be outliers
    get_unique_random_points!(v, np, np - p)
    
    # Add noise to some random points
    sgn = sign(randn())
    
    for k = 1:np
            
        y = model(x[k], θSol) + randn() * std

        noise = 0.0
            
        if k in v 
            y = model(x[k], θSol)
            noise = (1.0 + 2 * rand()) * out_times * std * sgn
        end
            
        data[k, 1] = x[k]
        data[k, 2] = y + noise
        data[k, 3] = noise

    end

    return data, θSol, v

end

"""

    generate_clustered_noisy_data(model::Function, n::Int, np::Int,
        p::Int, x_interval::Tuple{Float64,Float64},
        cluster_interval::Tuple{Float64, Float64}; kwargs...)

    generate_clustered_noisy_data(model::Function, n::Int,
        np::Int, p::Int, θSol::Vector{Float64},
        x_interval::Tuple{Float64,Float64},
        cluster_interval::Tuple{Float64, Float64}; kwargs...)

Generate a test set with clustered outliers.

The arguments and optional arguments are the same for
[`generate_noisy_data!`](@ref), with exception of tuple
`cluster_interval` which is the interval to generate the clustered
outliers.

It returns a tuple `(data, θSol, outliers)` where

  - `data`: (`np` x `3`) array, where each row contains `x` and
    `model(x, θSol)`. The same array given as argument
  - `θSol`: `n`-dimensional vector with the exact solution.
  - `outliers`: the outliers of this data set. The same vector given
    as argument.

"""
function generate_clustered_noisy_data(model::Function, n::Int,
    np::Int, p::Int, x_interval::Tuple{Float64,Float64},
    cluster_interval::Tuple{Float64, Float64}; kwargs...)

    data = Array{Float64, 2}(undef, np, 3)

    v = Vector{Int}(undef, np - p)

    return generate_clustered_noisy_data!(data, v, model, n, np, p,
               x_interval, cluster_interval; kwargs...)

end

generate_clustered_noisy_data(model::Function, n::Int, np::Int,
    p::Int, θSol::Vector{Float64}, x_interval::Tuple{Float64,Float64},
    cluster_interval::Tuple{Float64, Float64}; kwargs...) =
    generate_clustered_noisy_data(model, n, np, p, x_interval,
    cluster_interval, θSol=θSol; kwargs...)

"""

    generate_clustered_noisy_data!(data::Array{Float64, 2},
        v::Vector{Int}, model::Function, n::Int, np::Int, p::Int,
        x_interval::Tuple{Float64,Float64},
        cluster_interval::Tuple{Float64, Float64}; kwargs...)

Generate a test set with clustered outliers. This version overwrites
the content of (`np` x `3`) matrix `data` and vector `v` with integer
indices to the position of outliers in `data`.

The arguments and optional arguments are the same for
[`generate_noisy_data!`](@ref), with exception of tuple
`cluster_interval` which is the interval to generate the clustered
outliers.

It returns a tuple `(data, θSol, outliers)` where

  - `data`: (`np` x `3`) array, where each row contains `x` and
    `model(x, θSol)`. The same array given as argument
  - `θSol`: `n`-dimensional vector with the exact solution.
  - `outliers`: the outliers of this data set. The same vector given
    as argument.

"""
function generate_clustered_noisy_data!(data::Array{Float64, 2},
            v::Vector{Int}, model::Function, n::Int, np::Int,
            p::Int, x_interval::Tuple{Float64,Float64},
            cluster_interval::Tuple{Float64, Float64}; kwargs...)

    if (np - p > 0) &&
        !(x_interval[1] <= cluster_interval[1] <
          cluster_interval[2] <= x_interval[2])

        error("Bad interval for clustered data generation.")

    end
    
    interval_len = x_interval[2] - x_interval[1]
 
    fr1 = (cluster_interval[1] - x_interval[1]) / interval_len
    fr2 = (cluster_interval[2] - cluster_interval[1]) / interval_len
    fr3 = (x_interval[2] - cluster_interval[2]) / interval_len

    # Distribute data in a proportional way
    
    # Interval 2 will contain the clustered outliers
    np2 = max(np - p, Int(round(fr2 * np)))
    np1 = min(np - np2, Int(round(fr1 * np)))
    # Interval 3 will contain the remaining points
    np3 = max(0, np - np1 - np2)

    @debug("Clustered points: $(np1), $(np2), $(np3).")

    # Avoid repetition in the borders of the intervals
    
    δ_c = (cluster_interval[2] - cluster_interval[1]) / (np2 + 2)
    
    # Generate data

    tmpv = Vector{Int}()
    
    tmpd, θSol, tmpv = generate_noisy_data!(@view(data[1:np1, :]),
        tmpv, model, n, np1, np1; x_interval=(x_interval[1],
        cluster_interval[1]), kwargs...)

    generate_noisy_data!(@view(data[np1 + 1:np1 + np2, :]), v, model,
        n, np2, np2 - (np - p); x_interval=(cluster_interval[1] + δ_c,
        cluster_interval[2] - δ_c), θSol=θSol, kwargs...)

    # Update the outlier number with the correct number
    map!((x) -> x + np1, v, v)
    
    generate_noisy_data!(@view(data[np1 + np2 + 1:np, :]), tmpv,
        model, n, np3, np3; x_interval=(cluster_interval[2],
        x_interval[2]), θSol=θSol, kwargs...)

    return data, θSol, v
    
end

"""
    function generate_uniform_noisy_data!(data::AbstractArray{Float64, 2},
        model::Function, n::Int, np::Int, p::Int;
        x_interval::Tuple{Number, Number}=(-10.0, 10.0),
        θSol::Vector{Float64}=10.0 * randn(Float64, n),
        std::Number=2, thck::Int=2, funcsize::Int=200)

Modify figure `data` given by a real matrix so that it includes points
associated with `model` and random uniform noise. **Attention**: this
function only works with `1`-dimensional models.

The parameters are

  - `data`: the matrix representing a figure
  - `model`: real-valued model given by a function `model(x, θ)`
  - `n`: dimension of the parameters of the model
  - `np`: number of points to be generated
  - `p`: number of trusted points that will define the correct points
    in the model

The function also accepts the following optional arguments:

  - `x_interval`: tuple representing the interval for the `x` variable
  - `θSol`: vector with the 'exact' parameters of the solution
  - `std`: error that will be added to the simulated 'correct' points
  - `thck`: thickness of the point in the image
  - `funcsize`: size (in pixels) that the function will use in the image.

"""
function generate_uniform_noisy_data!(data::AbstractArray{Float64, 2},
    model::Function, n::Int, np::Int, p::Int;
    x_interval::Tuple{Number, Number}=(-10.0, 10.0),
    θSol::Vector{Float64}=10.0 * randn(Float64, n),
    std::Number=2, thck::Int=2, funcsize::Int=minimum(size(data)))

    @assert(x_interval[1] <= x_interval[2],
            "Invalid interval for random number generation.")

    # Generate (x_i) where x_interval[1] <= x_i <= x_interval[2] (data)
    # Fix the problem of large interval with 1 element.
    x = (np == 1) ? sum(x_interval) / 2.0 : LinRange(x_interval[1], x_interval[2], p)
    
    h, w = size(data)

    # Correct points are black
    y = map(x -> (model(x, θSol) + randn() * std), x)
    dx = minimum(x)
    dy = minimum(y)

    scale = funcsize / max(maximum(x) - dx, maximum(y) - dy)
    xshift = Int(round(rand() * (w - funcsize)))
    yshift = Int(round(rand() * (h - funcsize)))
    
    # Noise is gray
    for k = 1:np - p

        j = 1 + Int(round(rand() * w))
        i = 1 + Int(round(rand() * h))

        data[max(1, i - thck):min(i + thck, h), max(1, j - thck):min(j + thck, w)] .= 0.0

    end

    for (tx, ty) in zip(x, y)

        j = xshift + Int(round((- dx + tx) * scale))
        i = h - yshift - Int(round((- dy + ty) * scale))

        ((i < 1) || (i > h) || (j < 1) || (j > w)) && continue
        
        data[max(1, i - thck):min(i + thck, h), max(1, j - thck):min(j + thck, w)] .= 0.0
    end

    return data

end
