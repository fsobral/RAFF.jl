@testset "Parallel tests" begin

    # addprocs(2)

    # @everywhere using Logging
    # @everywhere using Base.CoreLogging

    # @everywhere global_logger(ConsoleLogger(stdout, Logging.Debug))

    
    # @everywhere gmodel!(x, t_, g) = begin
    #     g[1] = exp(t_ * x[2])
    #     g[2] = t_ * x[1] * exp(t_ * x[2]);
    # end

    # @everywhere model(x, t) = x[1] * exp(t * x[2])

    n = 2

    np = 10

    p = 7

    # data, xSol = RAFF.generateNoisyData(model, n, np, p)

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
    
end
