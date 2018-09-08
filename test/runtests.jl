using RAFF

using Test
using DelimitedFiles

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
    
    for line in eachline(datf)

        nLines += 1

        (parse(Float64, split(line)[3]) != 0.0) && (nNoise += 1)

    end

    @test (np - p) >= nNoise
    @test np == nLines        
    
end

@testset "Simple tests" begin

    model(x, t) = x[1] * exp(t * x[2])

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

    conv, x, iter, p = LMlovo(model, data, 2, 18)

    @test conv == 1
    @test x ≈ answer atol=1.0e-5
    @test p == 18
    
    conv, x, iter, p = raff(model, data, 2)
    
    @test conv == 1
    @test x ≈ answer atol=1.0e-5
    @test p == 18

    @test_throws AssertionError LMlovo(model, data, 0, 1)
    @test_throws AssertionError LMlovo(model, data, 2, -1)

    conv, x, iter, p = LMlovo(model, data, 2, 0)

    @test conv == 1
    @test iter == 1    
    
end

@testset "Generated test set" begin

    dir = "../examples/files/"

    # Iterate over a list of small problems and solutions
    for prob in eachline(dir * "list.dat")

        # Ignore blank lines
        (length(strip(prob)) == 0) && continue
        
        dname, sname = split(prob)

        # Data file
        data = readdlm(dir * dname)[:, [1, 2]]

        # Solution file
        fsol = open(dir * sname, "r")

        # Number of parameters
        n = Meta.parse(readline(fsol))

        # Solution vector
        answer = eval(Meta.parse(readline(fsol)))

        # Model function to fit data
        model = eval(Meta.parse(readline(fsol)))

        close(fsol)

        # Call raff
        conv, x, iter, p = raff(model, data, n)

        @test conv == 1
        @test x ≈ answer atol=1.0e-2

    end 

end
