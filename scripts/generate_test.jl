using RAFF
using Printf
# using ArgParse

# s = ArgParseSettings()
# @add_arg_table s begin

#     "np"
#     help = "Number of points"
#     arg_type = Int
#     required = true

#     "p"
#     help = "Number of trusted points"
#     arg_type = Int
#     required = true
    
#     "model"
#     help = "The model"
#     required = true

#     "-n"
#     help = "Number of parameters to find"
#     arg_type = Int
#     default = 2

#     "-s"
#     help = "Solution"
#     arg_type = Vector
    
# end

# Main program

# parsed_args = parse_args(ARGS, s)

# modelStr = parsed_args["model"]
    
# model = eval(Meta.parse(modelStr))

# println(parsed_args["s"])

model = (x, t) -> x[1] * t^3 + x[2] * t^2 + x[3] * t + x[4]

modelStr = "(x, t) -> x[1] * t^3 + x[2] * t^2 + x[3] * t + x[4]"

n = 2

np = 1000

p = 700

tMin = - 10.0

tMax =   10.0

xSol = [5.0, sqrt(7)]

fname = "L9"

# generateTestProblems("/tmp/prob.txt", "/tmp/sol.txt", model,
#                     modelStr, parsed_args["n"],
#                     parsed_args["np"], parsed_args["p"];
#                     )

try

    ffname = "/tmp/" * fname * ".txt"
    fsname = "/tmp/" * fname * "_sol.txt"
    
    generateTestProblems(ffname, fsname, model, modelStr, n, np, p;
                         tMin=tMin, tMax=tMax, xSol=xSol)

    @printf("Created problem and solution files.\n")

catch e

    @printf("Problems when creating test problem.\n")
    println(e)

end


