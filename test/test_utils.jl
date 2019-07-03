@testset "Util. tests" begin

    @testset "ftrusted" begin

        @test_throws MethodError RAFF.check_ftrusted(1, 15)
        
        @test RAFF.check_ftrusted(0.0, 15) == (0, 15)

        @test RAFF.check_ftrusted(0.3, 15) == (4, 15)

        @test RAFF.check_ftrusted(1.0, 15) == (15, 15)

        @test_throws ErrorException RAFF.check_ftrusted(-1.0, 15)

        @test RAFF.check_ftrusted((0.1, 0.4) , 15) == (2, 6)

        @test_throws ErrorException RAFF.check_ftrusted((0.4, 0.1) , 15)

        @test_throws ErrorException RAFF.check_ftrusted((0.1, 1.1) , 15)

        @test_throws ErrorException RAFF.check_ftrusted((-0.4, 0.1) , 15)
        
    end
    
    @testset "Sort function" begin

        ind = Vector{Int}(undef, 3)
        
        v = [3.0, 2.0, 1.0]

        pi, pv = RAFF.sort_fun!(v, ind, 2)

        @test(length(pi) == 2)
        @test(length(pv) == 2)
        @test(pi == [3, 2])
        @test(pv == [1.0, 2.0])
        @test(pi == ind[1:2])
        @test(pv == v[1:2])

        ind = Vector{Int}(undef, 4)
        
        v = [1.0, 3.0, 1.0, 1.0]

        pi, pv = RAFF.sort_fun!(v, ind, 3)

        @test(issubset(pi, [1, 3, 4]))
        @test(pv == ones(3))
        @test(pv == v[1:3])

        v = [1.0, 3.0, 1.0, 1.0]

        pi, pv = RAFF.sort_fun!(v, ind, 1)

        @test(pi == [1])
        @test(pv == [1.0])

        # Extreme behavior

        ind = zeros(Int, 3)
                          
        v = [3.0, 2.0, 1.0]

        pi, pv = RAFF.sort_fun!(v, ind, 0)

        @test(length(pi) == length(pv) == 0)
        @test(v == [3.0, 2.0, 1.0])
        @test(ind == zeros(3))

        pi, pv = RAFF.sort_fun!(v, ind, -1)

        @test(length(pi) == length(pv) == 0)
        @test(v == [3.0, 2.0, 1.0])
        @test(ind == zeros(3))
        
    end
    
end
