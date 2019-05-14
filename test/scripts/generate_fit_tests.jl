using RAFF
using Printf
using ArgParse

function mainf()

    s = ArgParseSettings()

    @add_arg_table s begin

        "np"
        help = "Number of points"
        arg_type = Int
        required = true

        "p"
        help = "Number of trusted points"
        arg_type = Int
        required = true

        "-m" "--model"
        help = "Model function"
        arg_type = String
        default = "linear"
        required = true

        "--fname"
        help = "File name"
        arg_type = String
        default = "output"
        
        "-s"
        help = "Solution"
        arg_type = Vector
        
    end

    # Main program

    parsed_args = parse_args(ARGS, s)

    n, model, modelStr = RAFF.model_list[parsed_args["m"]]

    xMin = - 0.0

    xMax =   30.0

    θSol = [1000.0, 5000.0, 0.2, 15.0]

    fname = parsed_args["fname"]

    try

        ffname = "/tmp/" * fname * ".txt"
        fsname = "/tmp/" * fname * "_sol.txt"
        
        generate_test_problems(ffname, fsname, model, modelStr, n,
                             parsed_args["np"], parsed_args["p"];
                             xMin=xMin, xMax=xMax, θSol=θSol, std=200.0)

        @printf("Created problem and solution files.\n")

    catch e

        @printf("Problems when creating test problem.\n")
        println(e)

    end

end

mainf()
