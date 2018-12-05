@testset "Multivariate" begin

    @testset "Simple test" begin

        data = [1.0 1.0   2.0
                0.0 0.0   4.0
                7.0 1.5  -4.5
                2.0 2.0 -17.0 # outlier
                0.0 8.6  -4.6]

        # This is the most complete way of defining the arguments for model and gmodel!
        model = (x::Vector{Float64}, t::Union{Vector{Float64}, SubArray}) -> x[1] *
            t[1] + x[2] * t[2] + x[3]

        gmodel! = (x::Vector{Float64}, t::Vector{Float64},
                   g::Union{Vector{Float64}, SubArray}) -> begin

            g[1] = t[1]
            g[2] = t[2]
            g[3] = 1.0
            
        end

        n = 3

        p = 4

        xsol = [- 1.0, - 1.0, 4.0]
        
        rout = lmlovo(model, gmodel!, data, n, p; ε=1.0e-8)

        @test rout.solution ≈ xsol atol=1.0e-4
        @test rout.outliers == [4]

        rout = raff(model, gmodel!, data, n)

        @test rout.solution ≈ xsol atol=1.0e-3
        @test rout.p == 4
        @test rout.outliers == [4]
        
        # Using automatic differentiation
        
        # This is the most complete way of defining the arguments for
        # model when automatic differentiation is being used
        model = (x, t::Union{Vector{Float64}, SubArray}) -> x[1] *
            t[1] + x[2] * t[2] + x[3]

        rout = lmlovo(model, data, n, p; ε=1.0e-8)

        @test rout.solution ≈ xsol atol=1.0e-4
        @test rout.outliers == [4]

        rout = raff(model, data, n)

        @test rout.solution ≈ xsol atol=1.0e-3
        @test rout.p == 4
        @test rout.outliers == [4]
        
    end
    
end
