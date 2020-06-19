using RAFF
using Printf
using ArgParse
using Random

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

        "--model"
        help = "Model function: linear, cubic, expon, logistic"
        arg_type = String
        default = "linear"

        "--fname"
        help = "File name"
        arg_type = String
        default = "output"
        
        "--sol"
        help = "Solution"
        nargs = '*'
        arg_type = Float64

        "--cint"
        help = "Cluster interval"
        nargs = 2
        arg_type = Float64

        "--i"
        help = "Generate i-th run."
        arg_type = Int
        default = 0

        "--std"
        help = "Standard deviation of the error."
        arg_type = Float64
        default = 200.0
        
    end

    # Main program

    parsed_args = parse_args(ARGS, s)

    n, model, modelStr = RAFF.model_list[parsed_args["model"]]

    x_interval = (- 0.0, 30.0)

    θSol = Vector{Float64}(undef, n)

    if length(parsed_args["sol"]) == 0

        randn!(θSol)

    elseif length(parsed_args["sol"]) == n

        θSol .= parsed_args["sol"]

    else

        error("Incorrect number of elements in solution")

    end

    fname = parsed_args["fname"]

    try

        ffname = "/tmp/" * fname * ".txt"
        fsname = "/tmp/" * fname * "_sol.txt"

        if parsed_args["i"] != 0

            Random.seed!(179424673 + parsed_args["i"])

        end

        if length(parsed_args["cint"]) == 2

            generate_test_problems(ffname, fsname, model, modelStr, n,
                parsed_args["np"], parsed_args["p"], x_interval,
                Tuple(parsed_args["cint"]); θSol=θSol, std=parsed_args["std"])

        else
            
            generate_test_problems(ffname, fsname, model, modelStr, n,
                parsed_args["np"], parsed_args["p"];
                x_interval=x_interval, θSol=θSol, std=parsed_args["std"])

        end
 
        @printf("Created problem and solution files.\n")

    catch e

        @printf("Problems when creating test problem.\n")
        println(e)

    end

end

mainf()
