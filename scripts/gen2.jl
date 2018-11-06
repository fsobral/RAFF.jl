using RAFF
using Printf
using ArgParse

const model_list = Dict(
    "linear" => (2, (x, t) -> x[1] * t + x[2],
                 "(x, t) -> x[1] * t + x[2]"),
    "cubic" => (2, (x, t) -> x[1] * t^3 + x[2] * t^2 + x[3] * t + x[4],
                "(x, t) -> x[1] * t^3 + x[2] * t^2 + x[3] * t + x[4]"),
    "expon" => (2, (x, t) -> x[1] * exp(- x[2] * t),
                "(x, t) -> x[1] * exp(- x[2] * t)")
)

function main()

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

    n, model, modelStr = model_list[parsed_args["m"]]

    tMin = - 0.0

    tMax =   30.0

    xSol = [2.0, 2.0]

    fname = parsed_args["fname"]

    try

        ffname = "/tmp/" * fname * ".txt"
        fsname = "/tmp/" * fname * "_sol.txt"
        
        generateTestProblems(ffname, fsname, model, modelStr, n,
                             parsed_args["np"], parsed_args["p"];
                             tMin=tMin, tMax=tMax, xSol=xSol)

        @printf("Created problem and solution files.\n")

    catch e

        @printf("Problems when creating test problem.\n")
        println(e)

    end

end

main()
