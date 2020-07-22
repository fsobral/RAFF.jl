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
        
    end

end
