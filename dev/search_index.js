var documenterSearchIndex = {"docs": [

{
    "location": "#",
    "page": "Overview",
    "title": "Overview",
    "category": "page",
    "text": ""
},

{
    "location": "#Overview-1",
    "page": "Overview",
    "title": "Overview",
    "category": "section",
    "text": "This page is devoted to document and make easier the use of RAFF- Robust  Algebraic Fitting Function. Our intent is provide a package  to determine fitting functions for a dataset with ability to detect possible  outliers of the dataset. All scripts were made in  Julia language, version 1.0. This package is not an implentation of classical least square solvers. In fact,  RAFF is a solver for Low Order Value Problem [1] which is a generalization of least square problems. "
},

{
    "location": "#Current-Status-1",
    "page": "Overview",
    "title": "Current Status",
    "category": "section",
    "text": "The current status of this project is beta quality, don\'t use for anything important.  We provide support to serial and parallel running. "
},

{
    "location": "#Developed-by-1",
    "page": "Overview",
    "title": "Developed by",
    "category": "section",
    "text": "This project was developed byFrancisco Sobral (Leader)\nEmerson Vitor Castelani\nRonaldo Lopes\nWesley Shirabayashi"
},

{
    "location": "#References-1",
    "page": "Overview",
    "title": "References",
    "category": "section",
    "text": "[1] Andreani, R., Martínez, J. M., Martínez, L., & Yano, F. S. (2009). Low order-value  optimization and applications. Journal of Global Optimization, 43(1), 1-22."
},

{
    "location": "tutorial/#",
    "page": "Tutorial",
    "title": "Tutorial",
    "category": "page",
    "text": ""
},

{
    "location": "tutorial/#Tutorial-1",
    "page": "Tutorial",
    "title": "Tutorial",
    "category": "section",
    "text": ""
},

{
    "location": "tutorial/#Installation-1",
    "page": "Tutorial",
    "title": "Installation",
    "category": "section",
    "text": "This package is supported just for Julia version 1.0. Consequently,  it uses package 3.0. Currently RAFF ins\'t in  Metadata.jl, so the  package can be installed with the Julia package manager. From the Julia REPL, type ] to enter into Pkg REPL mode and run:pkg> add https://github.com/fsobral/RAFF.jl#dev"
},

{
    "location": "tutorial/#Basic-usage-1",
    "page": "Tutorial",
    "title": "Basic usage",
    "category": "section",
    "text": "Just to ilustrate the potential and basic usage of RAFF, let us consider the following dataset given by an array A:A=left beginarraycc\n -20   50 \n -15   325\n -10   20 \n -05   125\n  00   10 \n  05   125\n  10   20 \n  15   325\n  20   50 \nendarrayrightLets suppose the first column of A as an experimental measure with  result given by second column. It is easy to see that the fitting  function in this case is accurate and given by phi(t)=t^2 +1Now lets perturb one result of second column of A. For example,  consider A_62=255. Assuming the model for fitting given byvarphi(xt)=x_1 t^2 +x_2 we have by classical least square as result x=[0.904329,1.3039]. But when we consider RAFF algorithm we obtain the correct answer x=[1.0,1.0]. Moreover, we have the number of possible outliers and which they are.In order to get run RAFF algorithm we need to setup using RAFFand define the dataset and model:A=[-2.0  5.0; \n  -1.5  3.25;\n  -1.0  2.0 ;\n  -0.5  1.25;\n   0.0  1.0 ;\n   0.5  2.55;\n   1.0  2.0 ;\n   1.5  3.25;\n   2.0  5.0 ;]\n\nmodel(x,t)=x[1]*t[1]^2+x[2]After that, we can run raff:raff(model,A,2)The number 2 above is the number of variables in model. The output is a RAFFoutput type, consequently is possible to handle with some atributes of this type. For example to acess only the parameters of solution, we can useoutput = raff(model,A,2)\noutput.solutionNote that RAFF algorithm detects and ignores possible outliers. In order to see which points are outiliers, we can acess the outliers atribute.  output.outliersMore details about RAFFoutput type and others options can be obtained in API section. By default RAFF uses automatic differentiation, more specifically ForwardDiff package. But is possible to define and handle RAFF with gradient vector of model. For example, considering the above example, we have, nabla varphi(xt)=t^21Programming this gradient and run raff we havegmodel(x,t,g)=begin\ng[1]=t[1]^2\ng[2]=1.0\nreturn g\nend\n\nraff(model,gmodel,A,2)"
},

{
    "location": "tutorial/#Multivariate-models-1",
    "page": "Tutorial",
    "title": "Multivariate models",
    "category": "section",
    "text": "RAFF supports fits to data sets of different sizes. To illustrate how this works, consider the following example:    data = [1.0 1.0   2.0\n            0.0 0.0   4.0\n            7.0 1.5  -4.5\n            2.0 2.0 -17.0 # outlier\n            0.0 8.6  -4.6]and the following model    model(x, t) = x[1]*t[1] + x[2] * t[2] + x[3]            Note that this model has two variables (t_1t_2). Naturally, this problem has one outlier (data[4,:]), so there are 4 trust points. Let\'s run RAFF and check the answer.     output = raff(model,data, 3)The right answer is [- 1.0, - 1.0, 4.0]. As we can note RAFF get a good fit for the data set.Handling the output follows the same pattern as the one-dimensional case. In order to get improvements in processing time, we can code the gradient vector of model to.    gmodel(x, t,g) = begin \n        g[1] = t[1]\n        g[2] = t[2]\n        g[3] = 1.0\n        return g\n    end    output = raff(model,data, 3)"
},

{
    "location": "tutorial/#Changing-some-options-1",
    "page": "Tutorial",
    "title": "Changing some options",
    "category": "section",
    "text": "Naturally, RAFF has options like precision of gradient stopping criteria and initial guess.     output = raff(model,data, 3;initguess=[0.5,0.5,0.5],ε=1.0e-4)RAFF is based on an optimization method. In this way, it is subject to stop at stationary points that are not global minimizers. For this reason, heuristics were implemented to find global minimizers. Such heuristics depend on random number generation. So if you want to run tests with more reliability this can be a useful strategy. To define in RAFF, say, 1000 different starting points, is enough to redefine the keyword argument MAXMS.    output = raff(model,data, 3;MAXMS=1,initguess=[0.5,0.5,0.5],ε=1.0e-10)"
},

{
    "location": "tutorial/#Parallel-running-1",
    "page": "Tutorial",
    "title": "Parallel running",
    "category": "section",
    "text": "RAFF can be run in a parallel or distributed environment, using the Distributed package and function praff. Let\'s use praff to solve the same problem from the beginning before. First, the Distributed has to be loaded and the number of workers is added. It is also possible to add the address of other machines.using Distributed\n\naddprocs(3) # Add 3 worker processesThis step can be replaced if Julia is initialized with the -p optionjulia -p 3Now we have to load RAFF and the fit function in all workers:@everywhere using RAFF\n\n@everywhere function gmodel!(x, t, g)\n    g[1] = t[1]^2\n    g[2] = 1.0\nend\n\n@everywhere function model(x, t)\n   x[1] * t[1]^2 + x[2]\nendthen, we call praff to solve the problem (note that we do not need to send the A matrix to all workers, since it will be sent by praff).A=[-2.0  5.0;\n  -1.5  3.25;\n  -1.0  2.0 ;\n  -0.5  1.25;\n   0.0  1.0 ;\n   0.5  2.55;\n   1.0  2.0 ;\n   1.5  3.25;\n   2.0  5.0 ;];\n\nn = 2\n\noutput = praff(model, gmodel!, A, n)\nRAFFOutput(1, [1.0, 0.999996], 6, 8, 4.0205772365906425e-11, [6])The true effectiveness of parallelism occurs when option MAXMS is set, which changes the number of random initial points that are tried for each subproblem solved. Better solutions can be achieved with higher values of MAXMSn = 2\n\noutput = praff(model, gmodel!, A, n; MAXMS=1000)\nRAFFOutput(1, [1.0, 1.0], 7, 8, 5.134133698545651e-13, [6])"
},

{
    "location": "api/#",
    "page": "API",
    "title": "API",
    "category": "page",
    "text": ""
},

{
    "location": "api/#Summary-1",
    "page": "API",
    "title": "Summary",
    "category": "section",
    "text": "There are three main RAFF structures: main functions: called by user; \nauxiliary functions: used like internal auxiliary function but can be modify user;\noutput type: type defined to manipulate output information."
},

{
    "location": "api/#RAFF.lmlovo",
    "page": "API",
    "title": "RAFF.lmlovo",
    "category": "function",
    "text": "lmlovo(model::Function [, x::Vector{Float64} = zeros(n)], data::Array{Float64, 2},\n       n::Int, p::Int [; kwargs...])\n\nlmlovo(model::Function, gmodel!::Function [, x::Vector{Float64} = zeros(n)],\n       data::Array{Float64,2}, n::Int, p::Int [; MAXITER::Int=200,\n       ε::Float64=10.0^-4])\n\nFit the n-parameter model model to the data given by matrix data. The strategy is based on the LOVO function, which means that only p (0 < p <= rows of data) points are trusted. The Levenberg-Marquardt algorithm is implemented in this version.\n\nMatriz data is the data to be fit. This matrix should be in the form\n\nt11 t12 ... t1N y1\nt21 t22 ... t2N y2\n:\n\nwhere N is the dimension of the argument of the model (i.e. dimension of t).\n\nIf \'x\' is provided, the it is used as the starting point.\n\nThe signature of function model should be given by\n\nmodel(x::Vector{Float64}, t::Union{Vector{Float64}, SubArray})\n\nwhere x is a n-dimensional vector of parameters and t is the argument. If the gradient of the model gmodel!\n\ngmodel!(x::Vector{Float64}, t::Union{Vector{Float64}, SubArray},\n        g::Vector{Float64})\n\nis not provided, then the function ForwardDiff.gradient! is called to compute it.  Note that this choice has an impact in the computational performance of the algorithm. In addition, if ForwardDiff is being used, then one MUST remove the signature of vector x from the model.\n\nThe optional arguments are\n\nMAXITER: maximum number of iterations\nε: tolerance for the gradient of the function\n\nReturns a RAFFOutput object.\n\n\n\n\n\n"
},

{
    "location": "api/#RAFF.raff",
    "page": "API",
    "title": "RAFF.raff",
    "category": "function",
    "text": "raff(model::Function, data::Array{Float64, 2}, n::Int; MAXMS::Int=1,\n     SEEDMS::Int=123456789, initguess=zeros(Float64, n))\n\nraff(model::Function, gmodel!::Function, data::Array{Float64, 2}, n::Int;\n     [MAXMS::Int=1, SEEDMS::Int=123456789, initguess=zeros(Float64, n),\n      kwargs...])\n\nRobust Algebric Fitting Function (RAFF) algorithm. This function uses a voting system to automatically find the number of trusted data points to fit the model.\n\nmodel: function to fit data. Its signature should be given by\nmodel(x, t)\nwhere x is a n-dimensional vector of parameters and t is the multidimensional argument\ngmodel!: gradient of the model function. Its signature should be given by\ngmodel!(x, t, g)\nwhere x is a n-dimensional vector of parameters, t is the multidimensional argument and the gradient is written in g.\ndata: data to be fit. This matrix should be in the form\nt11 t12 ... t1N y1\nt21 t22 ... t2N y2\n:\nwhere N is the dimension of the argument of the model (i.e. dimension of t).\nn: dimension of the parameter vector in the model function\n\nThe optional arguments are\n\nMAXMS: number of multistart points to be used\nSEEDMS: integer seed for random multistart points\ninitialguess: a good guess for the starting point and for generating random points in the multistart strategy\nε: gradient stopping criteria to lmlovo\n\nReturns a RAFFOutput object with the best parameter found.\n\n\n\n\n\n"
},

{
    "location": "api/#RAFF.praff",
    "page": "API",
    "title": "RAFF.praff",
    "category": "function",
    "text": "praff(model::Function, data::Array{Float64, 2}, n::Int; MAXMS::Int=1,\n      SEEDMS::Int=123456789, batches::Int=1, initguess=zeros(Float64, n),\n      ε=1.0e-4)\n\npraff(model::Function, gmodel!::Function, data::Array{Float64, 2}, n::Int;\n      MAXMS::Int=1, SEEDMS::Int=123456789, batches::Int=1,\n      initguess=zeros(Float64, n), ε::Float64=1.0e-4)\n\nMulticore distributed version of RAFF. See the description of the raff function for the main (non-optional) arguments. All the communication is performed by channels.\n\nThis function uses all available local workers to run RAFF algorithm. Note that this function does not use Tasks, so all the parallelism is based on the Distributed package.\n\nThe optional arguments are\n\nMAXMS: number of multistart points to be used\nSEEDMS: integer seed for random multistart points\nbatches: size of batches to be send to each worker\ninitguess: starting point to be used in the multistart procedure\nε: stopping tolerance\n\nReturns a RAFFOutput object containing the solution.\n\n\n\n\n\n"
},

{
    "location": "api/#RAFF.setRAFFOutputLevel",
    "page": "API",
    "title": "RAFF.setRAFFOutputLevel",
    "category": "function",
    "text": "setRAFFOutputLevel(level::LogLevel)\n\nSet the output level of raff and praff algorithms to the desired logging level. Options are (from highly verbose to just errors): Logging.Debug, Logging.Info, Logging.Warn and Logging.Error. The package Logging needs to be loaded.\n\nDefaults to Logging.Error.\n\n\n\n\n\n"
},

{
    "location": "api/#RAFF.setLMOutputLevel",
    "page": "API",
    "title": "RAFF.setLMOutputLevel",
    "category": "function",
    "text": "setLMOutputLevel(level::LogLevel)\n\nSet the output level of lmlovo algorithm to the desired logging level. Options are (from highly verbose to just errors): Logging.Debug, Logging.Info, Logging.Warn and Logging.Error. The package Logging needs to be loaded.\n\nDefaults to Logging.Error.\n\n\n\n\n\n"
},

{
    "location": "api/#Main-functions-1",
    "page": "API",
    "title": "Main functions",
    "category": "section",
    "text": "lmlovo\nraff\npraff\nsetRAFFOutputLevel\nsetLMOutputLevel"
},

{
    "location": "api/#RAFF.eliminate_local_min!",
    "page": "API",
    "title": "RAFF.eliminate_local_min!",
    "category": "function",
    "text": "eliminate_local_min!(sols::Vector{RAFFOutput})\n\nCheck if the function value of the solution found by smaller values of p is not greater when compared with larger ones. This certainly indicates that a local minimizer was found by the smaller p.\n\n\n\n\n\n"
},

{
    "location": "api/#RAFF.SortFun!",
    "page": "API",
    "title": "RAFF.SortFun!",
    "category": "function",
    "text": "This function is an auxiliary function. It finds the p smallest values of vector V and brings them to the first p positions. The indexes associated with the p smallest values are stored in ind.\n\n\n\n\n\n"
},

{
    "location": "api/#RAFF.update_best",
    "page": "API",
    "title": "RAFF.update_best",
    "category": "function",
    "text": "update_best(channel::RemoteChannel, bestx::SharedArray{Float64, 1})\n\nListen to a channel for results found by lmlovo. If there is an improvement for the objective function, the shared array bestx is updated.\n\nAttention: There might be an unstable state if there is a process   reading bestx while this function is updating it. This should not   be a problem, since it is used as a starting point.\n\nAttention 2: this function is currently out of use.\n\n\n\n\n\n"
},

{
    "location": "api/#RAFF.consume_tqueue",
    "page": "API",
    "title": "RAFF.consume_tqueue",
    "category": "function",
    "text": "function consume_tqueue(bqueue::RemoteChannel, tqueue::RemoteChannel,\n                        squeue::RemoteChannel, model::Function, gmodel!::Function,\n                        data::Array{Float64, 2}, n::Int, pliminf::Int,\n                        plimsup::Int, MAXMS::Int, seedMS::MersenneTwister)\n\nThis function represents one worker, which runs lmlovo in a multistart fashion.\n\nIt takes a job from the RemoteChannel tqueue and runs lmlovo function to it. It might run using a multistart strategy, if MAXMS>1. It sends the best results found for each value obtained in tqueue to channel squeue, which will be consumed by the main process. All the other arguments are the same for praff function.\n\n\n\n\n\n"
},

{
    "location": "api/#RAFF.check_and_close",
    "page": "API",
    "title": "RAFF.check_and_close",
    "category": "function",
    "text": "check_and_close(bqueue::RemoteChannel, tqueue::RemoteChannel,\n                squeue::RemoteChannel, futures::Vector{Future};\n                secs::Float64=0.1)\n\nCheck if there is at least one worker process in the vector of futures that has not prematurely finished. If there is no alive worker, close task, solution and best queues, tqueue, squeue and bqueue, respectively.\n\n\n\n\n\n"
},

{
    "location": "api/#RAFF.generateTestProblems",
    "page": "API",
    "title": "RAFF.generateTestProblems",
    "category": "function",
    "text": "generateTestProblems(datFilename::String, solFilename::String,\n                     model::Function, modelStr::String, n::Int,\n                     np::Int, p::Int)\n\nGenerate random data files for testing fitting problems.\n\ndatFilename and solFilename are strings with the name of the files for storing the random data and solution, respectively.\nmodel is the model function and modelStr is a string representing this model function, e.g.\n model = (x, t) -> x[1] * t[1] + x[2]\n modelStr = \"(x, t) -> x[1] * t[1] + x[2]\"\nwhere vector x represents the parameters (to be found) of the model and vector t are the variables of the model.\nn is the number of parameters\nnp is the number of points to be generated.\np is the number of trusted points to be used in the LOVO approach.\n\n\n\n\n\n"
},

{
    "location": "api/#RAFF.get_unique_random_points",
    "page": "API",
    "title": "RAFF.get_unique_random_points",
    "category": "function",
    "text": "get_unique_random_points(np::Int, npp::Int)\n\nChoose exactly npp unique random points from a set containing np points. This function is similar to rand(vector), but does not allow repetitions.\n\nReturn a vector with the selected points.\n\n\n\n\n\n"
},

{
    "location": "api/#RAFF.generateNoisyData",
    "page": "API",
    "title": "RAFF.generateNoisyData",
    "category": "function",
    "text": "generateNoisyData(model::Function, n::Int, np::Int, p::Int;\n                  tMin::Float64=-10.0, tMax::Float64=10.0,\n                  xSol::Vector{Float64}=10.0 * randn(Float64, n),\n                  std::Float64=200.0, outTimes::Float64=7.0)\n\ngenerateNoisyData(model::Function, n, np, p, tMin::Float64, tMax::Float64)\n\ngenerateNoisyData(model::Function, n::Int, np::Int, p::Int,\n                  xSol::Vector{Float64}, tMin::Float64, tMax::Float64)\n\nRandom generate a fitting one-dimensional data problem.\n\nThis function receives a model(x, t) function, the number of parameters n, the number of points np to be generated and the number of trusted points p. \n\nIf the n-dimensional vector xSol is provided, the the exact solution will not be random generated. The interval [tMin, tMax] for generating the values to evaluate model can also be provided.\n\nIt returns a tuple (data, xSol, outliers) where\n\ndata: (np x 2) array, where each row contains t and model(xSol, t).\nxSol: n-dimensional vector with the exact solution.\noutliers: the outliers of this data set\n\n\n\n\n\n"
},

{
    "location": "api/#Auxiliary-functions-1",
    "page": "API",
    "title": "Auxiliary functions",
    "category": "section",
    "text": "RAFF.eliminate_local_min!\nRAFF.SortFun!\nRAFF.update_best\nRAFF.consume_tqueue\nRAFF.check_and_close\nRAFF.generateTestProblems\nRAFF.get_unique_random_points\nRAFF.generateNoisyData"
},

{
    "location": "api/#RAFF.RAFFOutput",
    "page": "API",
    "title": "RAFF.RAFFOutput",
    "category": "type",
    "text": "This type defines the output file for the RAFF algorithm.\n\nRAFFOutput(status::Int, solution::Vector{Float64}, iter::Int,\n           p::Int, f::Float64, outliers::Vector{Int})\n\nwhere\n\nstatus: is 1 if converged and 0 if not\nsolution: vector with the parameters of the model\niter: number of iterations up to convergence\np: number of trusted points\nf: the residual value\noutliers: the possible outliers detected by the method, for the given p\n\nRAFFOutput()\n\nCreates a null version of output, equivalent to RAFFOutput(0, [], -1, 0, Inf, [])\n\nRAFFOuput(p::Int)\n\nCreates a null version of output for the given p.\n\n\n\n\n\n"
},

{
    "location": "api/#Output-type-1",
    "page": "API",
    "title": "Output type",
    "category": "section",
    "text": "RAFFOutput"
},

]}
