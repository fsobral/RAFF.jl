@testset "Util. tests" begin

    @testset "Sort function" begin

        ind = Vector{Int}(undef, 3)
        
        v = [3.0, 2.0, 1.0]

        pi, pv = RAFF.SortFun!(v, ind, 2)

        @test(length(pi) == 2)
        @test(length(pv) == 2)
        @test(pi == [3, 2])
        @test(pv == [1.0, 2.0])
        @test(pi == ind[1:2])
        @test(pv == v[1:2])

        ind = Vector{Int}(undef, 4)
        
        v = [1.0, 3.0, 1.0, 1.0]

        pi, pv = RAFF.SortFun!(v, ind, 3)

        @test(issubset(pi, [1, 3, 4]))
        @test(pv == ones(3))
        @test(pv == v[1:3])

        v = [1.0, 3.0, 1.0, 1.0]

        pi, pv = RAFF.SortFun!(v, ind, 1)

        @test(pi == [1])
        @test(pv == [1.0])

        # Extreme behavior

        ind = zeros(Int, 3)
                          
        v = [3.0, 2.0, 1.0]

        pi, pv = RAFF.SortFun!(v, ind, 0)

        @test(length(pi) == length(pv) == 0)
        @test(v == [3.0, 2.0, 1.0])
        @test(ind == zeros(3))

        pi, pv = RAFF.SortFun!(v, ind, -1)

        @test(length(pi) == length(pv) == 0)
        @test(v == [3.0, 2.0, 1.0])
        @test(ind == zeros(3))
        
    end
    
end

@testset "Random generator" begin

    modelStr = "(x, t) -> x[1] * t + x[2]"

    model = eval(Meta.parse(modelStr))

    n = 2

    np = 5

    p = 1

    datf = "test.dat"

    solf = "sol.dat"

    generateTestProblems(datf, solf, model, modelStr, n, np, p)

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

    end
    
end
