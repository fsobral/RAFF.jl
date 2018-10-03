@testset "Parallel tests" begin

    gmodel!(x, t_, g) = begin
        g[1] = 1.0
        g[2] = t_
    end

    model(x, t) = x[1]  + x[2] * t

    n = 2

    np = 10

    p = 7

    Random.seed!(123456789)
    
    data, xSol = RAFF.generateNoisyData(model, n, np, p)

    @debug(data)

    @testset "Updater" begin

        bqueue = RemoteChannel(() -> Channel{Vector{Float64}}(0))

        bestx = SharedArray{Float64, 1}(n)

        bestx .= 0.0

        fut = @async RAFF.update_best(bqueue, bestx)

        # Should update bestx

        newbest1 = ones(Float64, n)

        put!(bqueue, newbest1)

        sleep(0.5)

        @test bestx == newbest1

        # Should update bestx again

        newbest2 = rand(Float64, n)

        put!(bqueue, newbest2)

        sleep(0.5)

        @test bestx == (newbest1 + newbest2) / 2.0

        # Should not throw an error nor die with invalid element

        put!(bqueue, 10)

        @test !isready(bqueue)
        @test !istaskdone(fut)

        # Should finish

        close(bqueue)

        sleep(0.5)

        @test istaskdone(fut)
        
    end

    @testset "Consumer" begin

        bqueue = RemoteChannel(() -> Channel{Vector{Float64}}(0))

        tqueue = RemoteChannel(() -> Channel{Int}(0))

        bestx = SharedArray{Float64, 1}(n)

        bestx .= 0.0

        v = SharedArray{Float64, 2}(n, 6)
        vf = SharedArray{Float64, 1}(6)
        vs = SharedArray{Int, 1}(6)

        seedMS = MersenneTwister(1234)

        MAXMS = 1

        worker1 = @async RAFF.consume_tqueue(bqueue, tqueue,
                        bestx, v, vs, vf, model, gmodel!, data,
                        n, p - 2, np, MAXMS, seedMS)

        @test !istaskdone(worker1)

        # Should not do anything
        put!(tqueue, p - 3)

        # Should work, since the problem is easy
        put!(tqueue, p - 2)

        x = take!(bqueue)

        sleep(0.5)

        @test !istaskdone(worker1)
        @test x == v[:, 1]
        @test vs[1] == 1

        # Another test, with different p
        
        bestx .= xSol

        put!(tqueue, p)

        x2 = take!(bqueue)

        sleep(0.5)

        @test vs[3] == 1
        @test vf[3] ≈ 0.0 atol=1.0e-1
        @test v[:, 3] ≈ xSol atol=1.0e-1
        @test v[:, 3] == x2

        # Check if worker is alive when bqueue is closed

        close(bqueue)

        put!(tqueue, p)

        sleep(0.5)
        
        @test !istaskdone(worker1)

        # Check if worker finishes when tqueue is closed

        close(tqueue)

        sleep(0.1)

        @test istaskdone(worker1)

    end

    @testset "Consumer Multistart" begin

        bqueue = RemoteChannel(() -> Channel{Vector{Float64}}(4))

        tqueue = RemoteChannel(() -> Channel{Int}(0))

        bestx = SharedArray{Float64, 1}(n)

        bestx .= 0.0

        v = SharedArray{Float64, 2}(n, 6)
        vf = SharedArray{Float64, 1}(6)
        vs = SharedArray{Int, 1}(6)

        seedMS = MersenneTwister(1234)

        MAXMS = 1

        worker1 = @async RAFF.consume_tqueue(bqueue, tqueue,
                         bestx, v, vs, vf, model, gmodel!, data,
                         n, p - 2, np, MAXMS, seedMS)


        put!(tqueue, p)

        take!(bqueue)
        
        close(tqueue)

        sleep(0.1)

        @test istaskdone(worker1)

        # Save the objective function found and run multistart
        # strategy with the same initial random generator
        
        f1 = vf[3]

        tqueue = RemoteChannel(() -> Channel{Int}(0))

        MAXMS = 3

        seedMS = MersenneTwister(1234)

        worker = @async RAFF.consume_tqueue(bqueue, tqueue,
                        bestx, v, vs, vf, model, gmodel!, data,
                        n, p - 2, np, MAXMS, seedMS)

        put!(tqueue, p)

        close(tqueue)

        fetch(worker)

        close(bqueue)
        
        # Should find a better point
        
        @test f1 >= vf[3]
        @test istaskdone(worker)
        
    end
    
end
