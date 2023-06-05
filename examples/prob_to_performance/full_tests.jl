using DataFrames, DelimitedFiles, RAFF 
using LinearAlgebra, BenchmarkTools
using CSV


#######################################################
## Loading some useful things #########################
#######################################################

"""
    FitProbType

It is an immutable type used by main functions of this package

"""
struct FitProbType
    name::String
    data:: Array{Float64,2}
    npts:: Int
    nout:: Int
    model::Function
    dim:: Int
    cluster::Bool
    noise::Bool
    solution::Array{Float64,1}
    description::String
end

struct FitOutputType
    status::Bool
    solution::Vector{Float64}
    niter :: Int
    minimum :: Float64
    feval :: Int
    rval :: Int 
end

"""
    load_problem(filename::String)

This function is used to load a problem from a csv file and convert to FitProbType. It is an important function because FitProbType is the unique supported format in this package. 

# Examples
```
julia-repl
julia> load_problem("toy.csv")

returns a FitProbType
```
"""
function load_problem(filename::String)
    prob_matrix = readdlm(filename,':')
    return FitProbType(prob_matrix[1,2],eval(Meta.parse(prob_matrix[2,2])),prob_matrix[3,2],prob_matrix[4,2],eval(Meta.parse(prob_matrix[5,2])),prob_matrix[6,2],prob_matrix[7,2],prob_matrix[8,2],eval(Meta.parse(prob_matrix[9,2])),prob_matrix[10,2])
end


function solve(prob::FitProbType,θinit::Vector{Float64},method::String)
    m = prob.npts
    p = m-prob.nout
    A = prob.data
    model = prob.model
    dim = prob.dim


    if method == "LMlovo"
        res = RAFF.lmlovo(model, θinit, A, dim,p)
        return FitOutputType(res.status,res.solution,res.iter,res.f,res.nf,res.nj)
    end
    if method == "LMlovo_modified"
        res = RAFF.lmlovomodified(model, θinit, A, dim,p)
        return FitOutputType(res.status,res.solution,res.iter,res.f,res.nf,res.nj)
    end
    if method == "LMlovo_acelerated"
        res = RAFF.lmlovoacelerated(model, θinit, A, dim,p)
        return FitOutputType(res.status,res.solution,res.iter,res.f,res.nf,res.nj)
    end

    if method == "RAFF"
        res = RAFF.raff(model, A, dim;initguess = θinit)
        return FitOutputType(res.status,res.solution,res.iter,res.f,res.nf,res.nj)
    end

end

#######################################################
##### Performing ######################################
#######################################################

BenchmarkTools.DEFAULT_PARAMETERS.samples = 1000
BenchmarkTools.DEFAULT_PARAMETERS.seconds = 20
#function build_results(filename::String,method::String)
filename = "test_30_05.txt"
methods = ["LMlovo","LMlovo_modified"]
list = String.(readdlm(filename))
io = open("log_read_erros.txt", "w");
c = 1.0e-9
for methodname in methods
    # colocar em dataframe
    df = DataFrame(prob = String[],dim=Int64[],npts=Int64[],nout=Int64[], error= Float64[],minimum=Float64[],feval=Float64[],resval=Int64[],niter=Int64[],convergence=Bool[], minimum_PT = Float64[], maximum_PT = Float64[], mean_PT = Float64[], median_PT = Float64[], memory = Int64[],alloc = Int64[])
    #list = ["parabola_-0.5_-9.0_0.0_650_143.csv"]
    for file in list
        local prob = load_problem(file)
        println(" ====================================================================")
        println(" 🏁 Testing problem $(file) using $(methodname)")
        try  s = solve(prob,[1.0:1:prob.dim;],methodname)
            println(" 🏆 Problem free from read errors!")
            println(" 🕥 Performing benchmark! Please wait!")
            local b = @benchmark solve($(prob),[1.0:1:$(prob.dim);],$(methodname))
            display(b)
            push!(df,[file,prob.dim,prob.npts,prob.nout,norm(prob.solution-s.solution,Inf),s.minimum,s.feval,s.rval,s.niter,s.status,c*minimum(b.times),c*maximum(b.times),c*mean(b.times),c*median(b.times),b.memory,b.allocs])
        catch
            println(" 💣 Some read error in problem $(file) occurred")
            write(io," Error reading $(file) using $(methodname) occurred \n")
        end
       println(" ====================================================================")
    end
    CSV.write("$(methodname)_Output_Bench_Test.csv", df)
end
close(io)












#end
