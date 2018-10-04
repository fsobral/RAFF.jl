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

    @testset "Worker checker" begin

        nworkers = 3

        # Test if all workers are dead
        
        bqueue = RemoteChannel(() -> Channel{Vector{Float64}}(4))
        
        tqueue = RemoteChannel(() -> Channel{Int}(0))

        futures = Vector{Future}(undef, nworkers)

        for i = 1:nworkers

            futures[i] = @spawn error()

        end

        RAFF.check_and_close(bqueue, tqueue, futures)

        @test !isopen(bqueue)
        @test !isopen(tqueue)

        # Test if there is at least one live worker

        bqueue = RemoteChannel(() -> Channel{Vector{Float64}}(4))
        
        tqueue = RemoteChannel(() -> Channel{Int}(0))

        futures = Vector{Future}(undef, nworkers)

        for i = 1:nworkers - 1

            futures[i] = @spawn error()

        end

        futures[nworkers] = @spawn take!(tqueue)

        RAFF.check_and_close(bqueue, tqueue, futures)

        @test isopen(bqueue)
        @test isopen(tqueue)

        put!(tqueue, 1)

        sleep(0.1)

        close(bqueue)

        close(tqueue)

        @test fetch(futures[nworkers]) == 1
        @test !isopen(bqueue)
        @test !isopen(tqueue)
        
    end

    @testset "PRAFF" begin

        model(x, t) = x[1] * exp(t * x[2])
        
        gmodel!(x, t, g) = begin
            
            g[1] = exp(t * x[2])
            g[2] = t * x[1] * exp(t * x[2])

        end

        data = [-1.0   3.2974425414002564;
                -0.75  2.9099828292364025;
                -0.5    2.568050833375483;
                -0.25  2.2662969061336526;
                0.0                  2.0;
                0.25   1.764993805169191;
                0.5   1.5576015661428098;
                0.75  1.5745785575819442; #noise
                1.0   1.2130613194252668;
                1.25  1.0705228570379806;
                1.5   0.9447331054820294;
                1.75  0.8337240393570168;
                2.0   0.7357588823428847;
                2.25  0.6493049347166995;
                2.5   0.5730095937203802;
                2.75  0.5056791916094929;
                3.0  0.44626032029685964;
                3.25  0.5938233504083881; #noise 
                3.5   0.3475478869008902;
                3.75 0.30670993368985694;
                4.0   0.5706705664732254; #noise
                ]
        
        answer = [2.0, -0.5]

        # Regular test
        
        x, f, p = praff(model, data, 2)
        
        @test x ≈ answer atol=1.0e-5
        @test p == 18
        
        x, f, p = praff(model, gmodel!, data, 2)

        fgood = f
        
        @test x ≈ answer atol=1.0e-5
        @test p == 18

        # Multistart test

        x, f, p = praff(model, data, 2; MAXMS=2)
        
        @test x ≈ answer atol=1.0e-5
        @test p == 18
        @test f >= fgood
        
    end
    
end
