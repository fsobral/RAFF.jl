
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
    model(x,t) = Base.invokelatest(prob.model,x,t)
    dim = prob.dim


    if method == "lmlovo"
        res = RAFF.lmlovo(model, θinit, A, dim,p)
        return FitOutputType(res.status,res.solution,res.iter,res.f,res.nf,res.nj)
    end
    if method == "mlmlovo"
        res = RAFF.mlmlovo(model, θinit, A, dim,p)
        return FitOutputType(res.status,res.solution,res.iter,res.f,res.nf,res.nj)
    end

    if method == "RAFF"
        res = RAFF.raff(model, A, dim;initguess = θinit)
        return FitOutputType(res.status,res.solution,res.iter,res.f,res.nf,res.nj)
    end

end
function build_problem(probtype::String,limit::Vector{Float64},params::Vector{Float64})
    if probtype == "hiperplan"
        println("a limit vector is need to discretize the interval, for example, [-10.0,10.0]")
        println("params need to be setup as [coefs of A, b and c,size_prob,npts,nout]")
        nout = Int(params[end])
        npts = Int(params[end-1])
        size_prob = Int(params[end-2])
        dim = length(params[1:end-3])
        n = size_prob
        b = zeros(n)
        c = 0.0
        k = 1
        mb = "["
        for i=1:n-1
            mb = mb*" t[$(k)],"
            b[i] = params[k]
            k += 1
        end
        mb = mb*"t[$(k)]]"
        b[n] = params[k]
        k += 1
        mc = " t[$(k)]"
        c = params[k]
        h(x) = b'*x+c
        #display(A)
        model = "(x,t) -> $(mb)'*x + $(mc)"
        r = (limit[2]-limit[1])/(npts-1)
        x = rand(limit[1]:r:limit[2],npts)
        for i=1:n-1
            x = [x rand(limit[1]:r:limit[2],npts)]
        end
        #display(x)
        #global coefs = params[1:end-3]
        y = zeros(npts)
        #
        for i=1:npts
            y[i] = h(x[i,:])
        end

        k = 1
        iout = []
        while k<=nout
            i = rand([1:npts;])
            if i ∉ iout
                push!(iout,i)
                k = k+1
            end
        end

        for k = 1:nout
            #x[iout[k],:]=x[iout[k],:]+randn(size_prob)
            rd  = randn()
            y[iout[k]]=y[iout[k]]+sign(rd)+rd
        end

        
        FileMatrix = ["name :" "Hiperplan";"data :" [[x y]]; "npts :" npts;"nout :" nout; "model :" model;"dim :"  dim; "cluster :" "false"; "noise :" "false"; "solution :" [params[1:end-3]]; "description :" "none"]

        open("hiperplan_$(params[1])_$(params[end-3])_$(size_prob)_$(dim)_$(npts)_$(nout).csv", "w") do io
            writedlm(io, FileMatrix)
        end


    end



    if probtype == "quadratic"
        println("a limit vector is need to discretize the interval, for example, [-10.0,10.0]")
        println("params need to be setup as [coefs of A, b and c,size_prob,npts,nout]")
        nout = Int(params[end])
        npts = Int(params[end-1])
        size_prob = Int(params[end-2])
        dim = length(params[1:end-3])
        n = size_prob
        mA = "["
        sA = [" " for i=1:n,j=1:n]
        A = zeros(n,n)
        b = zeros(n)
        c = 0.0
        k = 1
        for i=1:n
            for j=i:n
                A[i,j] = params[k]
                A[j,i] = A[i,j]
                sA[i,j] = " t[$(k)] "
                sA[j,i] = sA[i,j]
                k += 1
            end
        end
        for i=1:n
            for j=1:n
                mA = mA*sA[i,j]
                if j==n
                    if i==n
                        mA = mA*"]"
                    else
                        mA = mA*";"
                    end
                end
            end
        end
        mb = "["
        for i=1:n-1
            mb = mb*" t[$(k)],"
            b[i] = params[k]
            k += 1
        end
        mb = mb*"t[$(k)]]"
        b[n] = params[k]
        k += 1
        mc = " t[$(k)]"
        c = params[k]
        m(x) = x'*A*x+b'*x+c
        #display(A)
        model = "(x,t) -> x'*$(mA)*x + $(mb)'*x + $(mc)"
        r = (limit[2]-limit[1])/(npts-1)
        x = rand(limit[1]:r:limit[2],npts)
        for i=1:n-1
            x = [x rand(limit[1]:r:limit[2],npts)]
        end
        #display(x)
        #global coefs = params[1:end-3]
        y = zeros(npts)
        #
        for i=1:npts
            y[i] = m(x[i,:])
        end

        k = 1
        iout = []
        while k<=nout
            i = rand([1:npts;])
            if i ∉ iout
                push!(iout,i)
                k = k+1
            end
        end

        for k = 1:nout
            #x[iout[k],:]=x[iout[k],:]+randn(size_prob)
            rd = randn()
            y[iout[k]]=y[iout[k]]+sign(rd)+rd
        end

        ## A(t) = begin
        #         mA = zeros(size_prob,size_prob)
        #         k = 1
        #         for i=1:size_prob
        #             for j=i:size_prob
        #                 mA[i,j] = t[k]
        #                 mA[j,i] = t[k]
        #                 k += 1
        #             end
        #         end
        #         return mA
        #     end
        # b(t) = begin
        #     vb = zeros(size_prob)
        #     k = Int(((size_prob^2+size_prob)/2)+1)
        #         for i=1:size_prob
        #             vb[i] = t[k]
        #             k += 1
        #         end
        #         return vb
        #     end
        #     c(t) = t[Int(((size_prob^2+size_prob)/2)+size_prob+1)]
        # 
        # model = "(x,t) -> x'*A(t)*x + b(t)'*x + c(t)"


        FileMatrix = ["name :" "Quadratic";"data :" [[x y]]; "npts :" npts;"nout :" nout; "model :" model;"dim :"  dim; "cluster :" "false"; "noise :" "false"; "solution :" [params[1:end-3]]; "description :" "none"]

        open("quadratic_$(params[1])_$(params[end-3])_$(size_prob)_$(dim)_$(npts)_$(nout).csv", "w") do io
            writedlm(io, FileMatrix)
        end


    end
    if probtype == "line"
        println("a limit vector is need to discretize the interval, for example, [-10.0,10.0]")
        println("params need to be setup as [a,b,npts,nout]")
        npts = Int(params[3])
        nout = Int(params[4])
        r = (limit[2]-limit[1])/(npts-1)
        x = [limit[1]:r:limit[2];]
        y = params[1]*x .+params[2]
        nout = Int(params[4])
        k = 1
        iout = []
        while k<=nout
            i = rand([1:npts;])
            if i ∉ iout
                push!(iout,i)
                k = k+1
            end
        end

        for k = 1:nout
            #x[iout[k]]=x[iout[k]]+randn()
            rd = randn()
            y[iout[k]]=y[iout[k]]+sign(rd)+rd
        end


        FileMatrix = ["name :" "Line";"data :" [[x y]]; "npts :" npts;"nout :" nout; "model :" "(x,t) -> t[1]*x[1] +t[2]";"dim :" 2; "cluster :" "false"; "noise :" "false"; "solution :" [params[1:2]]; "description :" "none"]

        open("line_$(params[1])_$(params[2])_$(npts)_$(nout).csv", "w") do io
            writedlm(io, FileMatrix)
        end
    end

    if probtype == "parabola"
        println("a limit vector is need to discretize the interval, for example, [-10.0,10.0]")
        println("params need to be setup as [a,b,c,npts,nout]")
        npts = Int(params[4])
        nout = Int(params[5])
        r = (limit[2]-limit[1])/(npts-1)
        x = [limit[1]:r:limit[2];]
        y = params[1]*x.^2 .+params[2]*x .+params[3]
        k = 1
        iout = []
        while k<=nout
            i = rand([1:npts;])
            if i ∉ iout
                push!(iout,i)
                k = k+1
            end
        end

        for k = 1:nout
            #x[iout[k]]=x[iout[k]]+randn()
            rd = randn()
            y[iout[k]]=y[iout[k]]+sign(rd)+rd
        end


        FileMatrix = ["name :" "Parabola";"data :" [[x y]]; "npts :" npts;"nout :" nout; "model :" "(x,t) -> t[1]*x[1]^2 +t[2]*x[1] + t[3]";"dim :" 3; "cluster :" "false"; "noise :" "false"; "solution :" [params[1:3]]; "description :" "none"]

        open("parabola_$(params[1])_$(params[2])_$(params[3])_$(npts)_$(nout).csv", "w") do io
            writedlm(io, FileMatrix)
        end
    end
    if probtype == "cubic"
        println("a limit vector is need to discretize the interval, for example, [-10.0,10.0]")
        println("params need to be setup as [a,b,c,d,npts,nout]")
        npts = Int(params[5])
        nout = Int(params[6])
        r = (limit[2]-limit[1])/(npts-1)
        x = [limit[1]:r:limit[2];]
        y = params[1]*x.^3 .+params[2]*x.^2 .+params[3]*x .+params[4]
        k = 1
        iout = []
        while k<=nout
            i = rand([1:npts;])
            if i ∉ iout
                push!(iout,i)
                k = k+1
            end
        end

        for k = 1:nout
            #x[iout[k]]=x[iout[k]]+randn()            
            rd = randn()
            y[iout[k]]=y[iout[k]]+sign(rd)+rd

        end


        FileMatrix = ["name :" "Cubic";"data :" [[x y]]; "npts :" npts;"nout :" nout; "model :" "(x,t) -> t[1]*x[1]^3 +t[2]*x[1]^2 + t[3]*x[1] +t[4]";"dim :" 4; "cluster :" "false"; "noise :" "false"; "solution :" [params[1:4]]; "description :" "none"]

        open("cubic_$(params[1])_$(params[2])_$(params[3])_$(params[4])_$(npts)_$(nout).csv", "w") do io
            writedlm(io, FileMatrix)
        end
    end

    if probtype == "sphere2D"
        println("params need to be setup as [center,radious,npts,nout]")
        c = [params[1],params[2]]
        r = params[3]
        npts = Int(params[4])
        x = zeros(npts)
        y = zeros(npts)
        θ = [0.0:2*π/(npts-1):2*π;]
        for k=1:npts
            x[k] = c[1]+r*cos(θ[k])
            y[k] = c[2]+r*sin(θ[k])
        end
        nout = Int(params[5])
        k = 1
        iout = []
        while k<=nout
            i = rand([1:npts;])
            if i ∉ iout
                push!(iout,i)
                k = k+1
            end
        end
        for k = 1:nout
            x[iout[k]]=x[iout[k]]+rand([0.25*r:0.1*(r);(1+0.25)*r])
            y[iout[k]]=y[iout[k]]+rand([0.25*r:0.1*(r);(1+0.25)*r])
        end
        FileMatrix = ["name :" "sphere2D";"data :" [[x y zeros(npts)]]; "npts :" npts;"nout :" nout; "model :" "(x,t) -> (x[1]-t[1])^2 + (x[2]-t[2])^2 - t[3]^2";"dim :" 3; "cluster :" "false"; "noise :" "false"; "solution :" [push!(c,r)]; "description :" "none"]

        open("sphere2D_$(c[1])_$(c[2])_$(c[3])_$(nout).csv", "w") do io
            writedlm(io, FileMatrix)
        end

    end
    if probtype == "gaussian"
        println("a limit vector is need to discretize the interval, for example, [-10.0,10.0]")
        println("params need to be setup as [a,b,c,npts,nout]")
        npts = Int(params[4])
        nout = Int(params[5])
        r = (limit[2]-limit[1])/(npts-1)
        x = [limit[1]:r:limit[2];]
        g(w) = params[1]*exp(-((w-params[2])^2)/(params[3]^2))
        y = g.(x)
        k = 1
        iout = []
        while k<=nout
            i = rand([1:npts;])
            if i ∉ iout
                push!(iout,i)
                k = k+1
            end
        end

        for k = 1:nout
            #x[iout[k]]=x[iout[k]]+randn()
            rd = randn()
            y[iout[k]]=y[iout[k]]+sign(rd)+rd
        end


        FileMatrix = ["name :" "gaussian";"data :" [[x y]]; "npts :" npts;"nout :" nout; "model :" "(x,t) -> t[1]*exp(-((x[1]-t[2])^2)/(t[3]^2))";"dim :" 3; "cluster :" "false"; "noise :" "false"; "solution :" [params[1:3]]; "description :" "none"]

        open("gaussian_$(params[1])_$(params[2])_$(params[3])_$(npts)_$(nout).csv", "w") do io
            writedlm(io, FileMatrix)
        end

    end 

    if probtype == "log"
        println("a limit vector is need to discretize the interval, for example, [-10.0,10.0]")
        println("params need to be setup as [a,b,c,npts,nout]")
        npts = Int(params[4])
        nout = Int(params[5])
        r = (limit[2]-limit[1])/(npts-1)
        x = [limit[1]:r:limit[2];]
        y = params[1]*log.((params[2]*x .+params[3]))
        nout = Int(params[5])
        k = 1
        iout = []
        while k<=nout
            i = rand([1:npts;])
            if i ∉ iout
                push!(iout,i)
                k = k+1
            end
        end

        for k = 1:nout
            #x[iout[k]]=x[iout[k]]+randn()
            y[iout[k]]=y[iout[k]]+randn()
        end


        FileMatrix = ["name :" "log";"data :" [[x y]]; "npts :" npts;"nout :" nout; "model :" "(x,t) -> t[1]*log((t[2]*x[1] +t[3]))";"dim :" 3; "cluster :" "false"; "noise :" "false"; "solution :" [params[1:3]]; "description :" "none"]

        open("log2_$(params[1])_$(params[2])_$(params[3])_$(npts)_$(nout).csv", "w") do io
            writedlm(io, FileMatrix)
        end
    end

    if probtype == "trig"
        println("a limit vector is need to discretize the interval, for example, [-10.0,10.0]")
        println("params need to be setup as [a,b,c,npts,nout]")
        npts = Int(params[4])
        nout = Int(params[5])
        r = (limit[2]-limit[1])/(npts-1)
        x = [limit[1]:r:limit[2];]
        y = params[1]*sin.(params[2]*x) .+ params[3]
        nout = Int(params[5])
        k = 1
        iout = []
        while k<=nout
            i = rand([1:npts;])
            if i ∉ iout
                push!(iout,i)
                k = k+1
            end
        end

        for k = 1:nout
           # x[iout[k]]=x[iout[k]]+randn()
            rd = randn()
            y[iout[k]]=y[iout[k]]+sign(rd)+rd
        end


        FileMatrix = ["name :" "trig";"data :" [[x y]]; "npts :" npts;"nout :" nout; "model :" "(x,t) -> t[1]*sin(t[2]*x[1]) +t[3]";"dim :" 3; "cluster :" "false"; "noise :" "false"; "solution :" [params[1:3]]; "description :" "none"]

        open("trig_$(params[1])_$(params[2])_$(params[3])_$(npts)_$(nout).csv", "w") do io
            writedlm(io, FileMatrix)
        end
    end


end



