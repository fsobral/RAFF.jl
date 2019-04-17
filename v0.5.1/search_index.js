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
    "text": "This page is devoted to document and make easier the use of RAFF- Robust Algebraic Fitting Function. Our intent is to provide a package to determine fitting functions for a dataset with ability to detect possible outliers of the dataset. All the code was made in Julia language, version 1.0.This package is not an implentation of classical least square solvers. In fact, RAFF is a solver for Low Order Value Problem [1] which is a generalization of least square problems."
},

{
    "location": "#Current-Status-1",
    "page": "Overview",
    "title": "Current Status",
    "category": "section",
    "text": "The current status of this project is beta quality, don\'t use for anything important.  We provide support to serial and parallel running."
},

{
    "location": "#Developed-by-1",
    "page": "Overview",
    "title": "Developed by",
    "category": "section",
    "text": "This project was developed by the optimization group at Department of Mathematics, State University of Maringá, BrazilFrancisco Sobral (Leader)\nEmerson Vitor Castelani\nRonaldo Lopes\nWesley Shirabayashi"
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
    "text": "This package is supported just for Julia version 1.0. Consequently,  it uses package 3.0. Currently RAFF is registered in General Julia Registers, so the  package can be installed using the Julia package manager. From the Julia REPL, type ] to enter into Pkg REPL mode and run:pkg> add RAFFIn what follows, we provide some simple examples on how to solve problems with RAFF. All the examples, and some other ones, are given in the examples/ directory as Julia scripts."
},

{
    "location": "tutorial/#Basic-usage-1",
    "page": "Tutorial",
    "title": "Basic usage",
    "category": "section",
    "text": "Just to illustrate the potential and basic usage of RAFF, let us consider the following data set given by an array A:A=left beginarraycc\n -20   50 \n -15   325\n -10   20 \n -05   125\n  00   10 \n  05   125\n  10   20 \n  15   325\n  20   50 \nendarrayrightLet\'s suppose the first column of A as an experimental measure with  result given by second column. It is easy to see that the fitting  function in this case is accurate and given by phi(t)=t^2 +1Now let\'s perturb one result of second column of A. For example,  consider A_62 = 255. Assuming the model for fitting given byvarphi(xt)=x_1 t^2 +x_2 we have by classical least squares as result x = [0.904329, 1.3039]. On the other hand, when we consider the RAFF algorithm we obtain the correct answer x=[1.0, 1.0]. Moreover, we also have a list of the possible outliers.In order to run RAFF algorithm we need to setup using RAFFand define the data set and model:A=[-2.0  5.0; \n   -1.5  3.25;\n   -1.0  2.0 ;\n   -0.5  1.25;\n    0.0  1.0 ;\n    0.5  2.55;\n    1.0  2.0 ;\n    1.5  3.25;\n    2.0  5.0 ;];\n\nmodel(x, t) = x[1] * t[1]^2 + x[2]After that, we can run method raff:raff(model, A, 2)The number 2 above is the number of variables in model, i. e., the number of parameters to adjust in the model. The output is a RAFFOutput type. For example to access only the parameters of solution, we can useoutput = raff(model, A, 2)\noutput.solutionNote that RAFF algorithm detects and ignores possible outliers. In order to see which points are outliers, we can access the outliers attribute.output.outliersMore details about RAFFOutput type and other options can be obtained in API section.By default RAFF uses automatic differentiation, more specifically ForwardDiff package. But is possible to call RAFF methods with gradient vector of model. For example, considering the above example, we have,nabla varphi(x t) = t^2 1Programming this gradient and run raff we havegmodel(x, t, g)=begin\n   g[1] = t[1]^2\n   g[2] = 1.0\n   return g\nend\n\nraff(model, gmodel, A, 2)Preliminary tests have shown that the use of explicit derivatives is 10 times faster than automatic differentiation."
},

{
    "location": "tutorial/#Multivariate-models-1",
    "page": "Tutorial",
    "title": "Multivariate models",
    "category": "section",
    "text": "RAFF supports the use of multivariate fitting functions to data sets of different dimensions. To illustrate how this works, consider the following example:data = [1.0 1.0    2.0\n        0.0 0.0    4.0\n        7.0 1.5   -4.5\n        2.0 2.0  -17.0 # outlier\n        0.0 8.6   -4.6]and the following modelmodel(x, t) = x[1] * t[1] + x[2] * t[2] + x[3]            Note that this model has two variables (t_1 t_2). Naturally, this problem has one outlier (data[4,:]), so there are 4 trust points. Let\'s run RAFF and check the answer.output = raff(model, data, 3)The right answer is [- 1.0, - 1.0, 4.0]. As we can note, RAFF get a good fit for the data set. Handling the output follows the same pattern as the one-dimensional case.In order to get improvements in processing time, we can code the gradient vector of model too:gmodel(x, t, g) = begin \n    g[1] = t[1]\n    g[2] = t[2]\n    g[3] = 1.0\n    return g\nendoutput = raff(model,data, 3)"
},

{
    "location": "tutorial/#Changing-some-options-1",
    "page": "Tutorial",
    "title": "Changing some options",
    "category": "section",
    "text": "Naturally, RAFF has options like precision of gradient stopping criteria and initial guess. output = raff(model,data, 3;initguess=[0.5,0.5,0.5],ε=1.0e-4)RAFF is based on an optimization method. In this way, it is subject to stopping at stationary points that are not global minimizers. For this reason, heuristics were implemented to find global minimizers. Such heuristics depend on random number generation. So, if you want to run tests with more reliability this can be a useful strategy. To define in RAFF, say, 1000 different starting points, is enough to redefine the keyword argument MAXMS.output = raff(model, data, 3; MAXMS=1, initguess=[0.5,0.5,0.5], ε=1.0e-10)In the above example, we have also changed the starting point for the method. Also, the stopping criterion was changed to 10^-10, which means high accuracy when solving the subproblems. See RAFF API for all the possible options that can be used."
},

{
    "location": "tutorial/#Parallel-running-1",
    "page": "Tutorial",
    "title": "Parallel running",
    "category": "section",
    "text": "RAFF can be run in a parallel or distributed environment, using the Distributed package and function praff. Let\'s use praff to solve the same problem from the beginning. First, the Distributed package has to be loaded and the number of workers has to be added. It is also possible to add the address of other machines.using Distributed\n\naddprocs(3) # Add 3 worker processesThis step can be replaced if Julia is initialized with the -p optionjulia -p 3Now we have to load RAFF and the fit function in all workers:@everywhere using RAFF\n\n@everywhere function model(x, t)\n   x[1] * t[1]^2 + x[2]\nend\n\n@everywhere function gmodel!(x, t, g)\n    g[1] = t[1]^2\n    g[2] = 1.0\nendthen, we call praff to solve the problem (note that we do not need to send the A matrix to all workers, since it will be automatically sent by praff).A=[-2.0  5.0;\n  -1.5  3.25;\n  -1.0  2.0 ;\n  -0.5  1.25;\n   0.0  1.0 ;\n   0.5  2.55;\n   1.0  2.0 ;\n   1.5  3.25;\n   2.0  5.0 ;];\n\nn = 2\n\noutput = praff(model, gmodel!, A, n)\nRAFFOutput(1, [1.0, 0.999996], 6, 8, 4.0205772365906425e-11, [6])The true effectiveness of parallelism occurs when option MAXMS is set, which changes the number of random initial points that are tried for each subproblem solved. Better solutions can be achieved with higher values of MAXMSn = 2\n\noutput = praff(model, gmodel!, A, n; MAXMS=1000)\nRAFFOutput(1, [1.0, 1.0], 7, 8, 5.134133698545651e-13, [6])"
},

{
    "location": "examples/#",
    "page": "Examples",
    "title": "Examples",
    "category": "page",
    "text": ""
},

{
    "location": "examples/#Examples-1",
    "page": "Examples",
    "title": "Examples",
    "category": "section",
    "text": "In addition to the examples given in the Tutorial, the examples/ directory contains another ways of using RAFF. Currently, we only provide an example on how to load a problem from file, solve it using RAFF and visually check the results.cubic.jl: this example solves a problem using a cubic model, with 4 parameters. The example also illustrates how to use RAFF.model_list utility structure in order to load pre-defined models."
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
    "text": "There are four main RAFF structures: Main functions: directly called by user; \nAuxiliary functions: used like internal auxiliary functions;\nRandom generation: used to generate random sets of data, in order to test RAFF\nOutput type: type defined to manipulate output information."
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
    "text": "raff(model::Function, data::Array{Float64, 2}, n::Int; MAXMS::Int=1,\n     SEEDMS::Int=123456789, initguess=zeros(Float64, n))\n\nraff(model::Function, gmodel!::Function, data::Array{Float64, 2}, n::Int;\n     [MAXMS::Int=1, SEEDMS::Int=123456789, initguess=zeros(Float64, n),\n      kwargs...])\n\nRobust Algebric Fitting Function (RAFF) algorithm. This function uses a voting system to automatically find the number of trusted data points to fit the model.\n\nmodel: function to fit data. Its signature should be given by\nmodel(x, t)\nwhere x is a n-dimensional vector of parameters and t is the multidimensional argument\ngmodel!: gradient of the model function. Its signature should be given by\ngmodel!(x, t, g)\nwhere x is a n-dimensional vector of parameters, t is the multidimensional argument and the gradient is written in g.\ndata: data to be fit. This matrix should be in the form\nt11 t12 ... t1N y1\nt21 t22 ... t2N y2\n:\nwhere N is the dimension of the argument of the model (i.e. dimension of t).\nn: dimension of the parameter vector in the model function\n\nThe optional arguments are\n\nMAXMS: number of multistart points to be used\nSEEDMS: integer seed for random multistart points\ninitialguess: a good guess for the starting point and for generating random points in the multistart strategy\nε: gradient stopping criteria to lmlovo\nnoutliers: integer describing the maximum expected number of outliers. The default is half.\n\nReturns a RAFFOutput object with the best parameter found.\n\n\n\n\n\n"
},

{
    "location": "api/#RAFF.praff",
    "page": "API",
    "title": "RAFF.praff",
    "category": "function",
    "text": "praff(model::Function, data::Array{Float64, 2}, n::Int; MAXMS::Int=1,\n      SEEDMS::Int=123456789, batches::Int=1, initguess=zeros(Float64, n),\n      ε=1.0e-4)\n\npraff(model::Function, gmodel!::Function, data::Array{Float64, 2}, n::Int;\n      MAXMS::Int=1, SEEDMS::Int=123456789, batches::Int=1,\n      initguess=zeros(Float64, n), ε::Float64=1.0e-4)\n\nMulticore distributed version of RAFF. See the description of the raff function for the main (non-optional) arguments. All the communication is performed by channels.\n\nThis function uses all available local workers to run RAFF algorithm. Note that this function does not use Tasks, so all the parallelism is based on the Distributed package.\n\nThe optional arguments are\n\nMAXMS: number of multistart points to be used\nSEEDMS: integer seed for random multistart points\nbatches: size of batches to be send to each worker\ninitguess: starting point to be used in the multistart procedure\nε: stopping tolerance\nnoutliers: integer describing the maximum expected number of outliers. The default is half.\n\nReturns a RAFFOutput object containing the solution.\n\n\n\n\n\n"
},

{
    "location": "api/#RAFF.set_raff_output_level",
    "page": "API",
    "title": "RAFF.set_raff_output_level",
    "category": "function",
    "text": "set_raff_output_level(level::LogLevel)\n\nSet the output level of raff and praff algorithms to the desired logging level. Options are (from highly verbose to just errors): Logging.Debug, Logging.Info, Logging.Warn and Logging.Error. The package Logging needs to be loaded.\n\nDefaults to Logging.Error.\n\n\n\n\n\n"
},

{
    "location": "api/#RAFF.set_lm_output_level",
    "page": "API",
    "title": "RAFF.set_lm_output_level",
    "category": "function",
    "text": "set_lm_output_level(level::LogLevel)\n\nSet the output level of lmlovo algorithm to the desired logging level. Options are (from highly verbose to just errors): Logging.Debug, Logging.Info, Logging.Warn and Logging.Error. The package Logging needs to be loaded.\n\nDefaults to Logging.Error.\n\n\n\n\n\n"
},

{
    "location": "api/#Main-functions-1",
    "page": "API",
    "title": "Main functions",
    "category": "section",
    "text": "lmlovo\nraff\npraff\nset_raff_output_level\nset_lm_output_level"
},

{
    "location": "api/#RAFF.eliminate_local_min!",
    "page": "API",
    "title": "RAFF.eliminate_local_min!",
    "category": "function",
    "text": "eliminate_local_min!(sols::Vector{RAFFOutput})\n\nCheck if the function value of the solution found by smaller values of p is not greater when compared with larger ones. This certainly indicates that a local minimizer was found by the smaller p.\n\n\n\n\n\n"
},

{
    "location": "api/#RAFF.sort_fun!",
    "page": "API",
    "title": "RAFF.sort_fun!",
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
    "location": "api/#Auxiliary-functions-1",
    "page": "API",
    "title": "Auxiliary functions",
    "category": "section",
    "text": "RAFF.eliminate_local_min!\nRAFF.sort_fun!\nRAFF.update_best\nRAFF.consume_tqueue\nRAFF.check_and_close"
},

{
    "location": "api/#RAFF.generate_test_problems",
    "page": "API",
    "title": "RAFF.generate_test_problems",
    "category": "function",
    "text": "generate_test_problems(datFilename::String, solFilename::String,\n                     model::Function, modelStr::String, n::Int,\n                     np::Int, p::Int)\n\nGenerate random data files for testing fitting problems.\n\ndatFilename and solFilename are strings with the name of the files for storing the random data and solution, respectively.\nmodel is the model function and modelStr is a string representing this model function, e.g.\n model = (x, t) -> x[1] * t[1] + x[2]\n modelStr = \"(x, t) -> x[1] * t[1] + x[2]\"\nwhere vector x represents the parameters (to be found) of the model and vector t are the variables of the model.\nn is the number of parameters\nnp is the number of points to be generated.\np is the number of trusted points to be used in the LOVO approach.\n\n\n\n\n\n"
},

{
    "location": "api/#RAFF.get_unique_random_points",
    "page": "API",
    "title": "RAFF.get_unique_random_points",
    "category": "function",
    "text": "get_unique_random_points(np::Int, npp::Int)\n\nChoose exactly npp unique random points from a set containing np points. This function is similar to rand(vector), but does not allow repetitions.\n\nReturn a vector with the selected points.\n\n\n\n\n\n"
},

{
    "location": "api/#RAFF.generate_noisy_data",
    "page": "API",
    "title": "RAFF.generate_noisy_data",
    "category": "function",
    "text": "generate_noisy_data(model::Function, n::Int, np::Int, p::Int;\n                  tMin::Float64=-10.0, tMax::Float64=10.0,\n                  xSol::Vector{Float64}=10.0 * randn(Float64, n),\n                  std::Float64=200.0, outTimes::Float64=7.0)\n\ngenerate_noisy_data(model::Function, n, np, p, tMin::Float64, tMax::Float64)\n\ngenerate_noisy_data(model::Function, n::Int, np::Int, p::Int,\n                  xSol::Vector{Float64}, tMin::Float64, tMax::Float64)\n\nRandom generate a fitting one-dimensional data problem.\n\nThis function receives a model(x, t) function, the number of parameters n, the number of points np to be generated and the number of trusted points p. \n\nIf the n-dimensional vector xSol is provided, the the exact solution will not be random generated. The interval [tMin, tMax] for generating the values to evaluate model can also be provided.\n\nIt returns a tuple (data, xSol, outliers) where\n\ndata: (np x 2) array, where each row contains t and model(xSol, t).\nxSol: n-dimensional vector with the exact solution.\noutliers: the outliers of this data set\n\n\n\n\n\n"
},

{
    "location": "api/#RAFF.model_list",
    "page": "API",
    "title": "RAFF.model_list",
    "category": "constant",
    "text": "This dictionary represents the list of models used in the generation of random tests. Return the tuple (n, model, model_str), where\n\nn is the number of parameters of the model\nmodel is the model of the form m(x, t), where x are the parameters and t are the variables\nmodel_str is the string representing the model, used to build random generated problems\n\n\n\n\n\n"
},

{
    "location": "api/#Random-generation-1",
    "page": "API",
    "title": "Random generation",
    "category": "section",
    "text": "RAFF.generate_test_problems\nRAFF.get_unique_random_points\nRAFF.generate_noisy_data\nRAFF.model_list"
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

{
    "location": "advanced/#",
    "page": "Advanced",
    "title": "Advanced",
    "category": "page",
    "text": ""
},

{
    "location": "advanced/#Advanced-usage-1",
    "page": "Advanced",
    "title": "Advanced usage",
    "category": "section",
    "text": ""
},

{
    "location": "advanced/#Test-Problems-1",
    "page": "Advanced",
    "title": "Test Problems",
    "category": "section",
    "text": "In order to develop a robust code, test problems regarding data fitting with outliers have been created. All the problems are available in directory test/test_problems. All the problems have the formbeginarrayccccc\n	k\n	t_11  t_12    t_1k  f_1\n	     \n	t_np1  t_np2    t_npk  f_np\nendarraywhere k is the number of variables of the model (not the number of parameters!), np is the number of data points to be adjusted, t_ij are the values selected for parameter j in experiment i and f_i is the result obtained for experiment i."
},

{
    "location": "advanced/#Script-files-1",
    "page": "Advanced",
    "title": "Script files",
    "category": "section",
    "text": "During the development and testing of RAFF several scripts and pieces of Julia code have been created. Those files are mostly related to the automated generation of test problems and visualization of the solutions. All those files are located in the test/scripts directory.We explain each file below, so maybe more advanced users can modify and re-use the code to generate their own problems.calc_ratio.jl: this script contains a function that generates several tests for the same model and solve each one with RAFF. Then it prints the ratio of outliers that have successfully been detected but the method.\ndraw.jl: This script draws the solution of a given problem and also its data, which was taken from a file. Examples of such files are located in test/test_problems.\nrun_raff.jl: this script simply loads a problem data from a file, selects a given model and runs RAFF. Examples of such files are located in test/test_problems.\nrun_lmlovo.jl: the same as before, but only for lmlovo function. Used mostly for testing.\ngenerate_fit_tests.jl: script for generating random test problem files, using the pre-defined models given by RAFF.model_list. This function cannot be called inside Julia, since it uses ArgParse package.\ngen_circle.jl: specif script for generating random test problems related to the detection of circles in the plane. It also provides functions to draw the problem and the solution, which differ from the draw.jl script above."
},

]}
