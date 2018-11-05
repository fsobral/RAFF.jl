@testset "Generated test set" begin

    dir = "../examples/files/"

    # Iterate over a list of small problems and solutions
    for prob in eachline(dir * "list.dat")

        # Ignore blank lines
        (length(strip(prob)) == 0) && continue
        
        dname, sname = split(prob)

        # Data file

        open(dir * dname, "r") do fp

            global N = parse(Int, readline(fp))
            
            global data = readdlm(fp)[:, [1, 2]]

        end

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
        x, f, p = raff(model, data, n)

        @test x ≈ answer atol=1.0e-2

    end 

end
