@testset "Solver $(solver)" for solver in [RAFF.lmlovo, RAFF.gnlslovo]

    model = (x, θ) -> θ[1] * x[1] + θ[2]
    gmodel! = (g, x, θ) -> begin
        g[1] = x[1]
        g[2] = 1.0
    end

    @testset "Simple cases ($(np), $p)" for (np, p) in [(10, 9), (100, 90)]
        
        data = Array{Float64, 2}(undef, np, 2)
        n = 2
        θsol = [-p, 0.5 * np]
        
        for (i, x) in enumerate(range(-10, stop=10, length=np))
            data[i, 1] = x
            data[i, 2] = model(x, θsol)
        end

        θ = [1.0, 1.0]

        data[1:np - p, 2] .+= 100
        
        r = solver(model, gmodel!, θ, data, n, p)

        @test(r.outliers == [1:np - p;])
        @test(r.solution ≈ θsol, atol=1.0e-4)
        @test(r.status == 1)
        @test(r.p == p)
        @test(r.nf >= r.iter)
        @test(r.nj >= r.iter)

        θ = [1.0, 1.0]

        r = solver(model, θ, data, n, p)

        @test(r.outliers == [1:np - p;])
        @test(r.solution ≈ θsol, atol=1.0e-4)
        @test(r.status == 1)
        @test(r.p == p)
        @test(r.nf >= r.iter)
        @test(r.nj >= r.iter)


    end

    @testset "Wrong arguments" begin

        np = 10
        p  = 9
        
        data = Array{Float64, 2}(undef, np, 2)
        n = 2
        θsol = [-p, 0.5 * np]
        
        for (i, x) in enumerate(range(-10, stop=10, length=np))
            data[i, 1] = x
            data[i, 2] = model(x, θsol)
        end

        θ = [1.0, 1.0]

        data[1:np - p, 2] .+= 100
        
        r = solver(model, gmodel!, θ, data, n, p)

        @test_throws(AssertionError, solver(model, θ, data, 0, 1))
        @test_throws(AssertionError, solver(model, θ, data, 2, -1))

    end

    # Test to check Issue #1
    
    @testset "Error in printing" begin

        m(x, θ) = θ[1] * x[1]^2 + θ[2]

        A = [ -2.0  5.00;
              -1.5  3.25;
              -1.0  2.00;
              -0.5  1.25;
              0.0  1.00;
              0.5  1.25;
              1.0  2.00;
              1.5  3.25;
              2.0  5.00 ]

        θ = [0.0, 0.0]

        # Changes log just for this test
        rout = with_logger(Logging.NullLogger()) do
            
            lmlovo(m, θ, A, 2, 4)

        end

        @test rout.status == 1
        @test rout.p == 4
        
    end

end
