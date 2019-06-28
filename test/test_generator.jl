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

        data, θSol1, v = RAFF.generate_noisy_data(model, n, np, p, θSol, -10.0, 10.0)

        @test all(data[:, 1] .>= -10.0)

        @test all(data[:, 1] .<=  10.0)

        @test θSol == θSol1

        # Test the memory efficient version

        data = Array{Float64}(undef, np, 3)

        v = Vector{Int64}(undef, np)
        
        data1, θSol1, v1 = RAFF.generate_noisy_data!(data, v, model, n, np, p)

        @test data == data1

        @test v == v1
        
    end

end
