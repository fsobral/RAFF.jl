@testset "Parallel tests" begin

    gmodel!(x, t_, g) = begin
        g[1] = 1.0
        g[2] = t_[1]
    end

    model(x, t) = x[1]  + x[2] * t[1]

    n = 2

    np = 10

    p = 7

    Random.seed!(123456789)
    
    data, xSol, = RAFF.generateNoisyData(model, n, np, p; std=0.0)

    # Remove outlier information
    data = data[:, 1:2]

    @testset "Updater" begin

        bqueue = RemoteChannel(() -> Channel{Vector{Float64}}(0))

        bestx = SharedArray{Float64, 1}(n)

        bestx .= 0.0

        fut = @async RAFF.update_best(bqueue, bestx)

        # Should update bestx

        newbest1 = ones(Float64, n)

        put!(bqueue, newbest1)

        sleep(0.1)

        @test bestx == newbest1

        # Should update bestx again

        newbest2 = rand(Float64, n)

        put!(bqueue, newbest2)

        sleep(0.1)

        @test bestx == (newbest1 + newbest2) / 2.0

        # Should not throw an error nor die with invalid element

        put!(bqueue, 10)

        @test !isready(bqueue)
        @test !istaskdone(fut)

        # Should finish

        close(bqueue)

        sleep(0.1)

        @test istaskdone(fut)
        
    end

    @testset "Consumer" begin

        bqueue = RemoteChannel(() -> Channel{Vector{Float64}}(0))

        tqueue = RemoteChannel(() -> Channel{UnitRange{Int}}(0))

        squeue = RemoteChannel(() -> Channel{RAFFOutput}(0))

        seedMS = MersenneTwister(1234)

        MAXMS = 1

        worker1 = @async @test begin
            RAFF.consume_tqueue(bqueue, tqueue, squeue,
                                model, gmodel!, data, n, p - 2, np,
                                MAXMS, seedMS)
            true
        end

        @test !istaskdone(worker1)

        # Should not do anything for an invalid interval
        put!(tqueue, p - 3:p)

        @test !isready(squeue)

        # Should work, since the problem is easy
        put!(tqueue, p - 2:p - 2)

        # Should return a vector with 1 solution pair
        rout = take!(squeue)

        @test !istaskdone(worker1)
        @test rout.p == p - 2
        @test rout.status == 1

        # Another test, with different p
        
        put!(tqueue, p:p)

        rout = take!(squeue)

        @test rout.p == p
        @test rout.status == 1
        @test rout.f ≈ 0.0 atol=1.0e-1
        @test rout.solution ≈ xSol atol=1.0e-1

        # Test with interval

        put!(tqueue, p - 1:p)

        rout = take!(squeue)

        @test rout.status == 1
        @test rout.p == p - 1

        rout = take!(squeue)

        @test rout.status == 1
        @test rout.p == p

        @test !isready(squeue)

        # Check if worker is alive when bqueue is closed

        close(bqueue)

        put!(tqueue, p:p)

        take!(squeue)

        @test !istaskdone(worker1)

        # Check if worker finishes when tqueue is closed

        close(tqueue)

        sleep(0.1)

        @test istaskdone(worker1)

        #  Worker should die if solution queue is closed

        tqueue = RemoteChannel(() -> Channel{UnitRange{Int}}(0))

        seedMS = MersenneTwister(1234)

        MAXMS = 1

        worker2 = @async @test begin
            RAFF.consume_tqueue(bqueue, tqueue, squeue,
                                model, gmodel!, data, n, p - 2, np,
                                MAXMS, seedMS)
            true
        end

        @test !istaskdone(worker2)

        close(squeue)

        put!(tqueue, p:p)

        sleep(0.1)

        @test !isready(tqueue)

        @test istaskdone(worker2)

    end

    @testset "Consumer Multistart" begin

        bqueue = RemoteChannel(() -> Channel{Vector{Float64}}(4))

        tqueue = RemoteChannel(() -> Channel{UnitRange{Int}}(0))

        squeue = RemoteChannel(() -> Channel{RAFFOutput}(0))

        seedMS = MersenneTwister(1234)

        MAXMS = 1

        # Check if the worker dies after closing the task queue

        worker1 = @async @test begin
            RAFF.consume_tqueue(bqueue, tqueue, squeue,
                                model, gmodel!, data, n, p - 2, np,
                                MAXMS, seedMS)
            true
        end

        put!(tqueue, p:p)

        close(tqueue)

        rout1 = take!(squeue)

        sleep(0.1)

        @test istaskdone(worker1)

        # Save the objective function found and run multistart
        # strategy with the same initial random generator
        
        tqueue = RemoteChannel(() -> Channel{UnitRange{Int}}(0))

        MAXMS = 3

        seedMS = MersenneTwister(1234)

        worker = @async @test begin
            RAFF.consume_tqueue(bqueue, tqueue, squeue,
                                model, gmodel!, data, n, p - 2, np,
                                MAXMS, seedMS)
            true
        end

        put!(tqueue, p:p)

        rout2 = take!(squeue)

        close(tqueue)

        fetch(worker)

        close(bqueue)
        
        # Should find a better point
        
        @test rout1.f >= rout2.f
        @test rout1.p == rout2.p
        @test istaskdone(worker)
        
    end

    @testset "Worker checker" begin

        nworkers = 3

        # Test if all workers are dead
        
        bqueue = RemoteChannel(() -> Channel{Vector{Float64}}(4))
        
        tqueue = RemoteChannel(() -> Channel{UnitRange{Int}}(0))

        squeue = RemoteChannel(() -> Channel{RAFFOutput}(0))

        futures = Vector{Future}(undef, nworkers)

        for i = 1:nworkers

            futures[i] = @spawn error()

        end

        RAFF.check_and_close(bqueue, tqueue, squeue, futures)

        @test !isopen(bqueue)
        @test !isopen(tqueue)
        @test !isopen(squeue)

        # Test if there is at least one live worker

        bqueue = RemoteChannel(() -> Channel{Vector{Float64}}(4))
        
        tqueue = RemoteChannel(() -> Channel{UnitRange{Int}}(0))

        squeue = RemoteChannel(() -> Channel{RAFFOutput}(0))

        futures = Vector{Future}(undef, nworkers)

        for i = 1:nworkers - 1

            futures[i] = @spawn error()

        end

        futures[nworkers] = @spawn take!(tqueue)

        RAFF.check_and_close(bqueue, tqueue, squeue, futures)

        @test isopen(bqueue)
        @test isopen(tqueue)
        @test isopen(squeue)

        put!(tqueue, 1:1)

        RAFF.check_and_close(bqueue, tqueue, squeue, futures)

        @test !isopen(bqueue)
        @test !isopen(tqueue)
        @test !isopen(squeue)

        # Ensure that this checker does not closes the queue if the
        # workers have finished their job very fast

        bqueue = RemoteChannel(() -> Channel{Vector{Float64}}(4))

        tqueue = RemoteChannel(() -> Channel{UnitRange{Int}}(0))

        squeue = RemoteChannel(() -> Channel{RAFFOutput}(0))

        futures = Vector{Future}(undef, nworkers)

        for i = 1:nworkers

            futures[i] = @spawn(nothing)

        end

        # Simulates the case where all the tasks have already been
        # taken
        close(tqueue)

        RAFF.check_and_close(bqueue, tqueue, squeue, futures)

        @test isopen(bqueue)
        @test isopen(squeue)
        
    end

    @testset "PRAFF" begin

        model(x, t) = x[1] * exp(t[1] * x[2])
        
        gmodel!(x, t, g) = begin
            
            g[1] = exp(t[1] * x[2])
            g[2] = t[1] * x[1] * exp(t[1] * x[2])

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
        
        rout = praff(model, data, 2)
        
        @test rout.solution ≈ answer atol=1.0e-5
        @test rout.p == 18
        
        rout = praff(model, gmodel!, data, 2)

        fgood = rout.f
        
        @test rout.solution ≈ answer atol=1.0e-5
        @test rout.p == 18

        # Multistart test

        rout = praff(model, data, 2; MAXMS=2)
        
        @test rout.solution ≈ answer atol=1.0e-5
        @test rout.p == 18
        @test rout.f <= fgood
        
    end
    
end
