using BenchmarkTools, RAFF, DelimitedFiles
#######################################################
##### Performing ######################################
#######################################################
function perform_raff()
    BenchmarkTools.DEFAULT_PARAMETERS.samples = 100
    BenchmarkTools.DEFAULT_PARAMETERS.seconds = 400
    #function build_results(filename::String,method::String)
    filename = "to_solve"
    c = 1.0e-9
    solver_list = [mlmlovo,lmlovo]
    problem_list = String.(readdlm(filename))
    prob = load_problem(problem_list[1])
    A = prob.data
    model(x, t) = Base.invokelatest(prob.model, x, t)
    dim = prob.dim
    raff(model, A, dim; inner_solver=solver_list[1])
    for solver in solver_list
        df = DataFrame(prob=String[], dim=Int64[], npts=Int64[], nout=Int64[], error=Float64[], minimum=Float64[], feval=Float64[], resval=Int64[], niter=Int64[], convergence=Bool[], minimum_PT=Float64[], maximum_PT=Float64[], mean_PT=Float64[], median_PT=Float64[], memory=Int64[], alloc=Int64[])
        for file ∈ problem_list
            prob = load_problem(file)
            A = prob.data
            model(x, t) = Base.invokelatest(prob.model, x, t)
            dim = prob.dim
            println(" 🏆 Problem free from read errors!")
            println(" 🏁 Testing problem $(file) using $(solver)")
            r = raff(model, A, dim; inner_solver=solver)
            s = FitOutputType(r.status,r.solution,r.iter,r.f,r.nf,r.nj)
            stio = open("status.txt", "a")
            if r.status == true
                println(" 🕥 Status = true, performing benchmark! Please wait!")
                b = @benchmark raff($(model), $(A), $(dim); inner_solver=$(solver))
                display(b)
                write(stio, "🥳 The Benchmark of problem $(file) was done! \n")
                close(stio)       
                push!(df, [file, prob.dim, prob.npts, prob.nout, norm(prob.solution - s.solution, Inf), s.minimum, s.feval, s.rval, s.niter, s.status, c * minimum(b.times), c * maximum(b.times), c * mean(b.times), c * median(b.times), b.memory, b.allocs])
            else
                println(" ✋ Status = false, the chosen method did not converge to this problem!")
                write(stio, "😢 The Benchmark of problem $(file) was not done! \n")
                close(stio)
                push!(df, [file, prob.dim, prob.npts, prob.nout, norm(prob.solution - s.solution, Inf), s.minimum, s.feval, s.rval, s.niter, s.status, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000])
            end
        end
        CSV.write("$(solver)_output_bench_test_in_$(filename).csv", df)
        #output = raff(model,A,dim;inner_solver = solver)
        #print("output from $(solver) \n")
        #display(output)
    end
end


function perform!()
    BenchmarkTools.DEFAULT_PARAMETERS.samples = 100
    BenchmarkTools.DEFAULT_PARAMETERS.seconds = 400
    #function build_results(filename::String,method::String)
    filename = "to_solve"
    methods = ["mlmlovo"]
    list = String.(readdlm(filename))
    #remain_list = copy(list)
    prob = load_problem(list[1])
    solve(prob, rand(prob.dim), methods[1])
    io = open("log_read_erros.txt", "w")
    c = 1.0e-9
    #Base.invokelatest(prob.model)
    for methodname in methods
        # colocar em dataframe
        df = DataFrame(prob=String[], dim=Int64[], npts=Int64[], nout=Int64[], error=Float64[], minimum=Float64[], feval=Float64[], resval=Int64[], niter=Int64[], convergence=Bool[], minimum_PT=Float64[], maximum_PT=Float64[], mean_PT=Float64[], median_PT=Float64[], memory=Int64[], alloc=Int64[])
        #list = ["parabola_-0.5_-9.0_0.0_650_143.csv"]
        for file in list
            try
                global prob = load_problem(file)
                println(" ====================================================================")
                println(" 🏆 Problem free from read errors!")
                println(" 🏁 Testing problem $(file) using $(methodname)")
                s = solve(prob, 0.5 .* [1.0:1:prob.dim;], methodname)
                stio = open("status.txt", "a")
                if s.status == true
                    println(" 🕥 Status = true, performing benchmark! Please wait!")
                    b = @benchmark solve($(prob), 0.5 .* [1.0:1:prob.dim;], $(methodname))
                    display(b)
                    write(stio, "🥳 The Benchmark of problem $(file) was done! \n")
                    close(stio)
                    push!(df, [file, prob.dim, prob.npts, prob.nout, norm(prob.solution - s.solution, Inf), s.minimum, s.feval, s.rval, s.niter, s.status, c * minimum(b.times), c * maximum(b.times), c * mean(b.times), c * median(b.times), b.memory, b.allocs])
                else
                    println(" ✋ Status = false, the chosen method did not converge to this problem!")
                    write(stio, "😢 The Benchmark of problem $(file) was not done! \n")
                    close(stio)
                    push!(df, [file, prob.dim, prob.npts, prob.nout, norm(prob.solution - s.solution, Inf), s.minimum, s.feval, s.rval, s.niter, s.status, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000])
                end
            catch

                println(" 💣 Some error in problem $(file) occurred")
                write(io, " Error reading $(file) or calling benchmark using $(methodname) occurred \n")
            end
            println(" ====================================================================")
        end
        CSV.write("$(methodname)_Output_Bench_Test_in_$(filename).csv", df)
    end
    close(io)
end

function perform_all(folder_list)
    #folder_list = ["parabola","gaussian","trigonometric","quadratic","cubic","hiperplan"]
    #folder_list = ["quadratic"]
    for flr in folder_list
        cd(flr)
        perform!()
        cd("..")
    end

end
