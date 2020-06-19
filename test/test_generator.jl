@testset "Random generator" begin

    modelStr = "(x, θ) -> θ[1] * x[1] + θ[2]"

    model = eval(Meta.parse(modelStr))

    n = 2

    np = 5

    p = 1

    datf = "test.dat"

    solf = "sol.dat"

    generate_test_problems(datf, solf, model, modelStr, n, np, p)

    open(solf, "r") do fp

        @test n == parse(Int, readline(fp))
        @test n == length(split(readline(fp)))
        @test modelStr == readline(fp)

    end

    nLines = 0
    nNoise = 0

    open(datf, "r") do fp

        # Check the dimension of the problem
        @test 1 == parse(Int, readline(fp))
        
        for line in eachline(fp)

            nLines += 1

            (parse(Int, split(line)[3]) != 0.0) && (nNoise += 1)

        end

    end

    @test (np - p) == nNoise
    @test np == nLines

    @testset "Unique" begin

        @test length(RAFF.get_unique_random_points(5,  2)) == 2

        @test length(RAFF.get_unique_random_points(0, -1)) == 0

        @test length(RAFF.get_unique_random_points(5,  0)) == 0

        @test length(RAFF.get_unique_random_points(1,  2)) == 1

        v = RAFF.get_unique_random_points(5, 5)

        @test length(v) == 5

        cnt = 0
        
        for i = 1:length(v)
            
            (i in v) && (cnt += 1)

        end

        @test cnt == 5

        @test length(RAFF.get_unique_random_points(1, -1)) == 0
        
        @test length(RAFF.get_unique_random_points(-1, 1)) == 0

        #

        @test_throws(DimensionMismatch,
                     RAFF.get_unique_random_points!(
                         Vector{Int}(undef, 1), 3, 2))

    end

    @testset "Noisy" begin

        n, model, model_s = RAFF.model_list["linear"]

        θSol = [1.0, 1.0]
        
        np = 5

        p = 4

        # No outliers
        
        data, θSol1, v = RAFF.generate_noisy_data(model, n, np, np)

        @test length(v) == 0

        @test length(θSol1) == n

        @test size(data) == (np, 3)

        @test all(data[:, 3] .== 0.0)

        # All outliers
        
        data, θSol1, v = RAFF.generate_noisy_data(model, n, np, 0)

        @test length(v) == np

        @test length(θSol1) == n

        @test size(data) == (np, 3)

        @test all(data[:, 3] .!= 0.0)

        #

        data, θSol1, v = RAFF.generate_noisy_data(model, n, np, p)

        @test length(v) == np - p

        @test length(θSol1) == n

        @test size(data) == (np, 3)

        @test sum(data[:, 3] .!= 0.0) == np - p

        # Test if the solution provided is maintained and also the interval

        xMin = -10.0

        xMax =  10.0
        
        data, θSol1, v = RAFF.generate_noisy_data(model, n, np, p, θSol, xMin, xMax)

        @test all(data[:, 1] .>= xMin)

        @test all(data[:, 1] .<= xMax)

        @test θSol == θSol1

        # Test if the solution provided is maintained and also the interval

        x_interval = (-10.0, 10.0)
        
        data, θSol1, v = RAFF.generate_noisy_data(model, n, np, p, θSol, x_interval)

        @test all(data[:, 1] .>= x_interval[1])

        @test all(data[:, 1] .<= x_interval[2])

        @test θSol == θSol1

        # Test the memory efficient version

        data = Array{Float64}(undef, np, 3)

        v = Vector{Int64}(undef, np)
        
        data1, θSol1, v1 = RAFF.generate_noisy_data!(data, v, model, n, np, p)

        @test data == data1

        @test v == v1
        
    end

    @testset "Cluster" begin

        n, model, model_s = RAFF.model_list["linear"]

        θSol = [1.0, 1.0]
        
        np = 10

        p = 7

        x_int = (-10.0, 10.0)

        c_int = (0.0, 5.0)

        data, θSol1, v = generate_clustered_noisy_data(model, n, np,
            p, x_int, c_int)

        @test size(data) == (np, 3)

        @test length(θSol1) == n

        @test length(v) == np - p

        out_ind = findall(abs.(data[:, 3]) .> 0.0)

        @test length(out_ind) == length(v)

        @test all(data[v, 1] .>= c_int[1])

        @test all(data[v, 1] .<= c_int[2])

        # This loop checks if the points are generated in order and
        # there is no repetition between groups of points
        is_ordered = true
        for i = 1:np - 1
            (data[i, 1] >= data[i + 1, 1]) && (is_ordered = false)
        end
        
        @test is_ordered

        # Cluster interval at the beginning

        x_int = (-10.0, 10.0)

        c_int = (-10.0, 0.0)

        data, θSol1, v = generate_clustered_noisy_data(model, n, np,
            p, x_int, c_int)

        @test length(v) == np - p
        
        @test any(data[:, 1] .> 0.0)
        
        # Non enclosing cluster interval
        
        x_int = (-10.0, 10.0)

        c_int = (-11.0, 0.0)

        @test_throws ErrorException generate_clustered_noisy_data(model, n, np,
            p, x_int, c_int)

        # Singleton cluster interval with np - p > 1

        np = 10

        p = 5
        
        x_int = (-10.0, 10.0)

        c_int = (1.0, 1.0)

        @test_throws ErrorException generate_clustered_noisy_data(model, n, np,
            p, x_int, c_int)

        # Singleton cluster interval with no outliers

        np = 10

        p = 10
        
        x_int = (-10.0, 10.0)
        
        c_int = (1.0, 1.0)

        data, θSol1, v = generate_clustered_noisy_data(model, n, np,
            p, x_int, c_int)

        @test length(v) == 0

        # Cluster with only one element

        np = 10

        p = 8

        x_int = (1.0, 30.0)

        c_int = (5.0, 10.0)

        data, θSol1, v = generate_clustered_noisy_data(model, n, np,
            p, x_int, c_int)

        @test length(findall(data[:, 1] .<= 5.0)) == 1
        
    end

    @testset "Rand Interv." begin

        n = 2
        
        x = Vector{Float64}(undef, n)

        l = [1.0, -5.0]
        u = [2.0, -1.0]
        
        interval = [i for i in zip(l, u)]
        
        RAFF.interval_rand!(x, interval)

        @test all(x .>= l)
        
        @test all(x .<= u)

        n = 5

        x = zeros(Float64, n)

        RAFF.interval_rand!(x, interval)

        @test all(x[1:2] .>= l)
        
        @test all(x[1:2] .<= u)

        @test all(x[3:n] .== 0.0)

        n = 1

        x = zeros(Float64, n)

        @test_throws ErrorException RAFF.interval_rand!(x, interval)

        # Bad interval
        
        n = 2
        
        x = Vector{Float64}(undef, n)

        l = [ 1.0, -5.0]
        u = [-1.0, -1.0]
        
        interval = [i for i in zip(l, u)]
        
        @test_throws ErrorException RAFF.interval_rand!(x, interval)

        
    end

    @testset "Model list $(type)" for (type, (n, model, model_str)) in RAFF.model_list

        # TODO: Maybe we need to get the dimension of the model?
        x = (type == "circle" || type == "ellipse") ? rand(2) : rand()

        θ = rand(n)

        model2 = eval(Meta.parse(model_str))

        @test model(x, θ) ≈ model2(x, θ)
        
    end

    @testset "Circle generator" begin

        data, v = generate_circle("circle_remove1.dat", 2, 1)

        @test(size(data) == (2, 4))
        @test(length(v) == 1)

        open("circle_remove1.dat", "r") do fp

            @test(parse(Int, readline(fp)) == 2)

            dataf = DelimitedFiles.readdlm(fp)

            @test(dataf ≈ data, atol=1.0e-15)
            @test(dataf[v[1], 4] == 1)

        end

        rm("./circle_remove1.dat")

        # Test bad interval

        @test_throws(ErrorException,
                     generate_circle("circle_remove2.dat", 5, 2;
                                     interval=[0.0, π]))
        @test_throws(SystemError, open("circle_remove2.dat"))

        # Test good interval

        data, v = generate_circle("circle_remove3.dat", 5, 4;
                     interval=[0.0, π/4, π / 2, 3 * π / 4, π ])

        @test(size(data) == (5, 4))
        @test(length(v) == 1)

        rm("circle_remove3.dat")

    end

    @testset "Circle uniform generator" begin

        data, v = generate_ncircle("circle_remove1.dat", 2, 1)

        @test(size(data) == (2, 4))
        @test(length(v) == 1)

        open("circle_remove1.dat", "r") do fp

            @test(parse(Int, readline(fp)) == 2)

            dataf = DelimitedFiles.readdlm(fp)

            @test(dataf ≈ data, atol=1.0e-14)
            @test(dataf[v[1], 4] > 0)

        end

        rm("./circle_remove1.dat")

        # Test bad interval

        @test_throws(ErrorException,
                     generate_circle("circle_remove2.dat", 5, 2;
                                     interval=[0.0, π]))
        @test_throws(SystemError, open("circle_remove2.dat"))

        # Test uniform distribution

        np = 100
        p  = 10
        θSol = [7.0, 7.0, 5.0]
        
        data, v = generate_circle("circle_remove3.dat", np, p;
                                  θSol=θSol)

        @test(size(data) == (100, 4))
        @test(length(v) == np - p)

        open("circle_remove3.dat") do fp

            @test(parse(Int, readline(fp)) == 2)

            dataf = DelimitedFiles.readdlm(fp)
            @test(dataf ≈ data, atol=1.0e-14)
            @test(all(dataf[v, 4] .> 0))

            @test(all(θSol[1] - 2 * θSol[3] .<= dataf[v, 1] .<= θSol[1] + 2 * θSol[3]))
            @test(all(θSol[2] - 2 * θSol[3] .<= dataf[v, 1] .<= θSol[2] + 2 * θSol[3]))

        end

        rm("circle_remove3.dat")

    end

end
