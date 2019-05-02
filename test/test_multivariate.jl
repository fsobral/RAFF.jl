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

    @testset "Circle" begin

        data = [2.0              	1.0             	0.0
                1.766044443118978	1.6427876096865393	-3.3306690738754696e-16
                1.1736481776669305	1.9848077530122081	2.220446049250313e-16
                0.5000000000000002	1.8660254037844388	0.0
                0.06030737921409168	1.3420201433256689	-1.1102230246251565e-16
                0.06030737921409157	0.6579798566743313	0.0
                0.49999999999999956	0.13397459621556151	0.0
                3.17364817766693	0.015192246987791869	0.0 # noise
                1.766044443118978	0.3572123903134604	0.0]

        xsol = [1.0, 1.0, 1.0]

        model = (x, t::Union{Vector{Float64}, SubArray}) -> (t[1] - x[1])^2 + (t[2] - x[2])^2 - x[3]

        n = 3

        p = 8
        
        rout = lmlovo(model, [0.0, 0.0, 2.0], data, n, p; ε=1.0e-8)

        @test rout.solution ≈ xsol atol=1.0e-4
        @test rout.outliers == [8]

        rout = raff(model, data, n; ε=1.0e-8)

        @test rout.p == 8
        @test rout.outliers == [8]
        @test rout.solution ≈ xsol atol=1.0e-4

    end
    
end
