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
    "text": "This page is devoted to document and make easier the use of RAFF- Robust Algebraic Fitting Function. Our intent is to provide a package to determine fitting functions for a dataset with ability to detect possible outliers of the dataset. All the code was made in Julia language, version 1.0.This package is not an implentation of classical least squares solvers. It is an optimization-based package, based on algorithms for Lower Order-Value Optimization (LOVO) which were introduced in [1] and revisited in [2] to fit the user-provided models to experimental data. Recently, a good review can be found in [3]. To find possible outliers, LOVO methods depend on the number of outliers as input information. RAFF differs in this point and has no dependence on this number of outliers to perform the fitting process. In order to find a robust adjustment, a voting system is used, which is also responsible for the detection of possible outliers."
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
    "text": "This project was developed by the optimization group at Department of Mathematics, State University of Maringá, Brazil.Francisco Sobral (Leader)\nEmerson Vitor Castelani\nRonaldo Lopes\nWesley ShirabayashiThe authors of this package were sponsored by Fundação Araucária, project number 002/17 - 47223."
},

{
    "location": "#References-1",
    "page": "Overview",
    "title": "References",
    "category": "section",
    "text": "[1] Andreani, R., Dunder, C. & Martínez, J.M. Math Meth Oper Res (2005) 61: 365. https://doi.org/10.1007/s001860400410[2] Andreani, R., Martínez, J.M., Martínez, L. et al. J Glob Optim (2009) 43: 1. https://doi.org/10.1007/s10898-008-9280-3[3] Martínez, J.M. TOP (2012) 20: 75. https://doi.org/10.1007/s11750-010-0169-1"
},

{
    "location": "#Citing-this-package-1",
    "page": "Overview",
    "title": "Citing this package",
    "category": "section",
    "text": "If you would like to cite this package, please useCastelani, E. V., Lopes, R., Shirabayashi, W., & Sobral, F. N. C. (2019). RAFF.jl: Robust Algebraic Fitting Function in Julia. Journal of Open Source Software, 4(39), 1385. https://doi.org/10.21105/joss.01385BibTex"
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
    "text": "Just to illustrate the potential and basic usage of RAFF, let us consider the following data set given by an array A:A=left beginarraycc\n -20   50 \n -15   325\n -10   20 \n -05   125\n  00   10 \n  05   125\n  10   20 \n  15   325\n  20   50 \nendarrayrightLet\'s suppose the first column of A as an experimental measure with  result given by second column. It is easy to see that the fitting  function in this case is accurate and given by phi(x) = x^2 + 1Now let\'s perturb one result of second column of A. For example,  consider A_62 = 255. Assuming the model for fitting given byvarphi(x theta) = theta_1 x^2 + theta_2 we have by classical least squares as result θ = [0.904329, 1.3039]. On the other hand, when we consider RAFF algorithm we obtain the correct answer θ = [1.0, 1.0]. Moreover, we also have a list of the possible outliers.In order to run RAFF algorithm we need to setup using RAFFand define the data set and model:A=[-2.0  5.0; \n   -1.5  3.25;\n   -1.0  2.0 ;\n   -0.5  1.25;\n    0.0  1.0 ;\n    0.5  2.55;\n    1.0  2.0 ;\n    1.5  3.25;\n    2.0  5.0 ;];\n\nmodel(x, θ) = θ[1] * x[1]^2 + θ[2]After that, we can run method raff:raff(model, A, 2)The number 2 above is the number of variables in model, i. e., the number of parameters to adjust in the model. The output is a RAFFOutput type. For example, to access only the parameters of solution, we can useoutput = raff(model, A, 2)\noutput.solutionNote that RAFF algorithm detects and ignores possible outliers. In order to see which points are outliers, we can access the outliers attribute.output.outliersMore details about RAFFOutput type and other options can be obtained in API section.By default RAFF uses automatic differentiation, more specifically ForwardDiff.jl package. But is possible to call RAFF methods with gradient vector of model. For example, considering the above example, we have,nabla varphi(x theta) = x^2 1Programming this gradient and running RAFF we havegmodel!(g, x, θ) = begin\n   g[1] = x[1]^2\n   g[2] = 1.0\nend\n\nraff(model, gmodel!, A, 2)Preliminary tests have shown that the use of explicit derivatives is 10 times faster than automatic differentiation."
},

{
    "location": "tutorial/#Multivariate-models-1",
    "page": "Tutorial",
    "title": "Multivariate models",
    "category": "section",
    "text": "RAFF supports the use of multivariate fitting functions to data sets of different dimensions. To illustrate how this works, consider the following example:data = [1.0 1.0    2.0\n        0.0 0.0    4.0\n        7.0 1.5   -4.5\n        2.0 2.0  -17.0 # outlier\n        0.0 8.6   -4.6]and the following modelmodel(x, θ) = θ[1] * x[1] + θ[2] * x[2] + θ[3]Note that this model has two variables (x_1 x_2) and three parameters (theta_1 theta_2 theta_3). This problem has one outlier (data[4,:]), so there are 4 trusted points. Let\'s run RAFF and check the answer.output = raff(model, data, 3)The right answer is [-1.0, -1.0, 4.0]. As we can note, RAFF get a good fit for the data set. Handling the output follows the same pattern as the one-dimensional case.In order to get improvements in processing time, we can code the gradient vector of model too:gmodel!(g, x, θ) = begin \n    g[1] = x[1]\n    g[2] = x[2]\n    g[3] = 1.0\nendoutput = raff(model, gmodel!, data, 3)"
},

{
    "location": "tutorial/#Changing-some-options-1",
    "page": "Tutorial",
    "title": "Changing some options",
    "category": "section",
    "text": "RAFF has tunning options like precision of gradient stopping criteria and initial guess.output = raff(model, data, 3; initguess=[0.5,0.5,0.5], ε=1.0e-4)RAFF is based on an optimization method. In this way, it is subject to stopping at stationary points that are not global minimizers. For this reason, heuristics were implemented to find global minimizers. Such heuristics depend on random number generation. So, if you want to run tests with more reliability this can be a useful strategy. To define in RAFF, say, 1000 different starting points, is enough to redefine the keyword argument MAXMS.output = raff(model, data, 3; MAXMS=1, initguess=[0.5,0.5,0.5], ε=1.0e-10)In the above example, we have also changed the starting point for the method. Also, the stopping criterion was changed to 10^-10, which means high accuracy when solving the subproblems. See RAFF API for all the possible options that can be used."
},

{
    "location": "tutorial/#Parallel-running-1",
    "page": "Tutorial",
    "title": "Parallel running",
    "category": "section",
    "text": "RAFF can be run in a parallel or distributed environment, using the Distributed package and function praff. Let\'s use praff to solve the same problem from the beginning. First, the Distributed package has to be loaded and the number of workers has to be added. It is also possible to add the address of other machines.using Distributed\n\naddprocs(3) # Add 3 worker processesThis step can be replaced if Julia is initialized with the -p optionjulia -p 3Now we have to load RAFF and the fit function in all workers:@everywhere using RAFF\n\n@everywhere function model(x, θ)\n   θ[1] * x[1]^2 + θ[2]\nend\n\n@everywhere function gmodel!(g, x, θ)\n    g[1] = x[1]^2\n    g[2] = 1.0\nendthen, we call praff to solve the problem (note that we do not need to send the A matrix to all workers, since it will be automatically sent by praff).A=[-2.0  5.0;\n  -1.5  3.25;\n  -1.0  2.0 ;\n  -0.5  1.25;\n   0.0  1.0 ;\n   0.5  2.55;\n   1.0  2.0 ;\n   1.5  3.25;\n   2.0  5.0 ;];\n\nn = 2\n\noutput = praff(model, gmodel!, A, n)\nRAFFOutput(1, [1.0, 0.999996], 6, 8, 4.0205772365906425e-11, [6])The true effectiveness of parallelism occurs when option MAXMS is set, which changes the number of random initial points that are tried for each subproblem solved. Better solutions can be achieved with higher values of MAXMSn = 2\n\noutput = praff(model, gmodel!, A, n; MAXMS=1000)\nRAFFOutput(1, [1.0, 1.0], 7, 8, 5.134133698545651e-13, [6])"
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
    "text": "In addition to the examples given in the Tutorial, the examples/ directory contains another ways of using RAFF. Currently, we only provide an example on how to load a problem from file, solve it using RAFF and visually check the results.cubic.jl: this example solves a problem using a cubic model, with 4 parameters. The example also illustrates how to use RAFF.model_list utility structure in order to load pre-defined models.draw_and_detect.jl: this nice example uses GtkReactive.jl to show a graphic application of RAFF.jl to the detection of circles drawn by the user. The user can also see the difference between the LOVO approach and the traditional least squares technique."
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
    "text": "lmlovo(model::Function [, θ::Vector{Float64} = zeros(n)], data::Array{Float64, 2},\n       n::Int, p::Int [; kwargs...])\n\nlmlovo(model::Function, gmodel!::Function [, θ::Vector{Float64} = zeros(n)],\n       data::Array{Float64,2}, n::Int, p::Int [; MAXITER::Int=200,\n       ε::Float64=10.0^-4])\n\nFit the n-parameter model model to the data given by matrix data. The strategy is based on the LOVO function, which means that only p (0 < p <= rows of data) points are trusted. The Levenberg-Marquardt algorithm is implemented in this version.\n\nMatriz data is the data to be fit. This matrix should be in the form\n\nx11 x12 ... x1N y1\nx21 x22 ... x2N y2\n:\n\nwhere N is the dimension of the argument of the model (i.e. dimension of x).\n\nIf θ is provided, then it is used as the starting point.\n\nThe signature of function model should be given by\n\nmodel(x::Union{Vector{Float64}, SubArray}, θ::Vector{Float64})\n\nwhere x are the variables and θ is a n-dimensional vector of parameters. If the gradient of the model gmodel!\n\ngmodel! = (g::SubArray, x::Union{Vector{Float64}, SubArray},\n           θ::Vector{Float64})\n\nis not provided, then the function ForwardDiff.gradient! is called to compute it.  Note that this choice has an impact in the computational performance of the algorithm. In addition, if ForwardDiff.jl is being used, then one MUST remove the signature of vector θ from function model.\n\nThe optional arguments are\n\nMAXITER: maximum number of iterations\nε: tolerance for the gradient of the function\n\nReturns a RAFFOutput object.\n\n\n\n\n\n"
},

{
    "location": "api/#RAFF.gnlslovo",
    "page": "API",
    "title": "RAFF.gnlslovo",
    "category": "function",
    "text": "gnlslovo(model, gmodel!, θ, data::Array{T, 2}, n, p;\n         ε::Number=1.0e-4, MAXITER=400, αls=2.0, dinc=2.0,\n         MAXLSITER=100) where {T<:Float64}\n\ngnlslovo(model, θ::Vector{Float64}, data::Array{Float64,2},\n         n::Int, p::Int; kwargs...)\n\ngnlslovo(model, gmodel!, data::Array{Float64,2}, n::Int,\n         p::Int; kwargs...)\n\ngnlslovo(model, data::Array{Float64,2}, n::Int, p::Int; kwargs...)\n\nLOVO Gauss-Newton with line-search described in\n\nR. Andreani, G. Cesar, R. M. Cesar-Jr., J. M. Martínez, and P. J. S. Silva, “Efficient curve detection using a {Gauss-Newton} method with applications in agriculture,” in Proc. 1st International Workshop on Computer Vision Applications for Developing Regions in Conjunction with ICCV 2007-CVDR-ICCV07, 2007.\n\nFit the n-parameter model model to the data given by matrix data. The strategy is based on the LOVO function, which means that only p (0 < p <= rows of data) points are trusted.\n\nMatriz data is the data to be fit. This matrix should be in the form\n\nx11 x12 ... x1N y1\nx21 x22 ... x2N y2\n:\n\nwhere N is the dimension of the argument of the model (i.e. dimension of x).\n\nIf θ is provided, then it is used as the starting point.\n\nThe signature of function model should be given by\n\nmodel(x::Union{Vector{Float64}, SubArray}, θ::Vector{Float64})\n\nwhere x are the variables and θ is a n-dimensional vector of parameters. If the gradient of the model gmodel!\n\ngmodel! = (g::SubArray, x::Union{Vector{Float64}, SubArray},\n           θ::Vector{Float64})\n\nis not provided, then the function ForwardDiff.gradient! is called to compute it.  Note that this choice has an impact in the computational performance of the algorithm. In addition, if ForwardDiff.jl is being used, then one MUST remove the signature of vector θ from function model.\n\nThe optional arguments are\n\nMAXITER: maximum number of iterations\nε: tolerance for the gradient of the function\nαls: number >1 to increase/decrease the parameter t in line-search\ndinc: number >1 to increase the diagonal of the J^T J matrix in order to escape from singularity\nMAXLSITER: maximum number of Linear System increases in diagonal before exiting. Also defines the maximum number of Line Search trials to satisfy Armijo (but does not exit in such case)\n\nReturns a RAFFOutput object.\n\n\n\n\n\n"
},

{
    "location": "api/#RAFF.raff",
    "page": "API",
    "title": "RAFF.raff",
    "category": "function",
    "text": "raff(model::Function, data::Array{Float64, 2}, n::Int; kwargs...)\n\nraff(model::Function, gmodel!::Function, data::Array{Float64, 2},\n    n::Int; MAXMS::Int=1, SEEDMS::Int=123456789,\n    initguess::Vector{Float64}=zeros(Float64, n),\n    noutliers::Int=-1, ftrusted::Union{Float64,\n    Tuple{Float64, Float64}}=0.5,\n    inner_solver::Function=lmlovo, inner_solver_params...)\n\nRobust Algebric Fitting Function (RAFF) algorithm. This function uses a voting system to automatically find the number of trusted data points to fit the model.\n\nmodel: function to fit data. Its signature should be given by\nmodel(x, θ)\nwhere x is the multidimensional argument and θ is the n-dimensional vector of parameters\ngmodel!: gradient of the model function. Its signature should be given by\ngmodel!(g, x, θ)\nwhere x is the multidimensional argument, θ is the n-dimensional vector of parameters and the gradient is written in g.\ndata: data to be fit. This matrix should be in the form\nx11 x12 ... x1N y1\nx21 x22 ... x2N y2\n:\nwhere N is the dimension of the argument of the model (i.e. dimension of x).\nn: dimension of the parameter vector in the model function\n\nThe optional arguments are\n\nMAXMS: number of multistart points to be used\nSEEDMS: integer seed for random multistart points\ninitialguess: a good guess for the starting point and for generating random points in the multistart strategy\nnoutliers: integer describing the maximum expected number of outliers. The default is half. Deprecated.\nftrusted: float describing the minimum expected percentage of trusted points. The default is half (0.5). Can also be a Tuple of the form (fmin, fmax) percentages of trusted points.\ninner_solver: solver to be used for the least square problems. By default, uses lmlovo. This function has the following mandatory parameters\ninner_solver(model, gmodel!, θ, data, n, p;\n             inner_solver_params...) = RAFFOutput\ninner_solver_params...: the remaining parameters will be sent as optional arguments to the inner_solver\n\nReturns a RAFFOutput object with the best parameter found.\n\n\n\n\n\n"
},

{
    "location": "api/#RAFF.praff",
    "page": "API",
    "title": "RAFF.praff",
    "category": "function",
    "text": "praff(model::Function, data::Array{Float64, 2}, n::Int; kwargs...)\n\npraff(model::Function, gmodel!::Function, data::Array{Float64, 2},\n    n::Int; MAXMS::Int=1, SEEDMS::Int=123456789, batches::Int=1,\n    initguess::Vector{Float64}=zeros(Float64, n),\n    noutliers::Int=-1, ftrusted::Union{Float64,\n    Tuple{Float64, Float64}}=0.5,\n    inner_solver::Function=lmlovo, inner_solver_params...)\n\nMulticore distributed version of RAFF. See the description of the raff function for the main (non-optional) arguments. All the communication is performed by channels.\n\nThis function uses all available local workers to run RAFF algorithm. Note that this function does not use Tasks, so all the parallelism is based on the Distributed package.\n\nThe optional arguments are\n\nMAXMS: number of multistart points to be used\nSEEDMS: integer seed for random multistart points\nbatches: size of batches to be send to each worker\ninitguess: starting point to be used in the multistart procedure\nnoutliers: integer describing the maximum expected number of outliers. The default is half. Deprecated.\nftrusted: float describing the minimum expected percentage of trusted points. The default is half (0.5). Can also be a Tuple of the form (fmin, fmax) percentages of trusted points.\ninner_solver: solver to be used for the least square problems. By default, uses lmlovo. This function has the following mandatory parameters\ninner_solver(model, gmodel!, θ, data, n, p;\n             inner_solver_params...) = RAFFOutput\ninner_solver_params...: the remaining parameters will be sent as optional arguments to the inner_solver\n\nReturns a RAFFOutput object containing the solution.\n\n\n\n\n\n"
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
    "text": "lmlovo\ngnlslovo\nraff\npraff\nset_raff_output_level\nset_lm_output_level"
},

{
    "location": "api/#Auxiliary-functions-1",
    "page": "API",
    "title": "Auxiliary functions",
    "category": "section",
    "text": "RAFF.voting_strategy\nRAFF.eliminate_local_min!\nRAFF.sort_fun!\nRAFF.update_best\nRAFF.consume_tqueue\nRAFF.check_and_close\nRAFF.check_ftrusted\nRAFF.interval_rand!"
},

{
    "location": "api/#RAFF.generate_test_problems",
    "page": "API",
    "title": "RAFF.generate_test_problems",
    "category": "function",
    "text": "generate_test_problems(datFilename::String, solFilename::String,\n    model::Function, modelStr::String, n::Int, np::Int, p::Int;\n    x_interval::Tuple{Float64, Float64}=(-10.0, 10.0),\n    θSol::Vector{Float64}=10.0 * randn(n), std::Float64=200.0,\n    out_times::Float64=7.0)\n\ngenerate_test_problems(datFilename::String, solFilename::String,\n    model::Function, modelStr::String, n::Int, np::Int, p::Int,\n    cluster_interval::Tuple{Float64, Float64};\n    x_interval::Tuple{Float64, Float64}=(-10.0, 10.0),\n    θSol::Vector{Float64}=10.0 * randn(n), std::Float64=200.0,\n    out_times::Float64=7.0)\n\nGenerate random data files for testing fitting problems.\n\ndatFilename and solFilename are strings with the name of the files for storing the random data and solution, respectively.\nmodel is the model function and modelStr is a string representing this model function, e.g.\n model = (x, θ) -> θ[1] * x[1] + θ[2]\n modelStr = \"(x, θ) -> θ[1] * x[1] + θ[2]\"\nwhere vector θ represents the parameters (to be found) of the model and vector x are the variables of the model.\nn is the number of parameters\nnp is the number of points to be generated.\np is the number of trusted points to be used in the LOVO approach.\n\nIf cluster_interval is provided, then generates outliers only in this interval.\n\nAdditional parameters:\n\nxMin, xMax: interval for generating points in one dimensional tests Deprecated\nx_interval: interval for generating points in one dimensional tests\nθSol: true solution, used for generating perturbed points\nstd: standard deviation\nout_times: deviation for outliers will be out_times * std.\n\n\n\n\n\n"
},

{
    "location": "api/#RAFF.get_unique_random_points",
    "page": "API",
    "title": "RAFF.get_unique_random_points",
    "category": "function",
    "text": "get_unique_random_points(np::Int, npp::Int)\n\nChoose exactly npp unique random points from a set containing np points. This function is similar to rand(vector), but does not allow repetitions.\n\nIf npp < np, returns all the np points. Note that this function is not very memory efficient, since the process of selecting unique elements involves creating several temporary vectors.\n\nReturn a vector with the selected points.\n\n\n\n\n\n"
},

{
    "location": "api/#RAFF.get_unique_random_points!",
    "page": "API",
    "title": "RAFF.get_unique_random_points!",
    "category": "function",
    "text": "get_unique_random_points!(v::Vector{Int}, np::Int, npp::Int)\n\nChoose exactly npp unique random points from a set containing np points. This function is similar to rand(vector), but does not allow repetitions.\n\nIf npp < np, returns all the np points. Note that this function is not very memory efficient, since the process of selecting unique elements involves creating several temporary vectors.\n\nReturn the vector v provided as argument filled with the selected points.\n\n\n\n\n\n"
},

{
    "location": "api/#RAFF.generate_noisy_data!",
    "page": "API",
    "title": "RAFF.generate_noisy_data!",
    "category": "function",
    "text": "generate_noisy_data!(data::AbstractArray{Float64, 2},\n    v::Vector{Int}, model::Function, n::Int, np::Int, p::Int;\n    x_interval::Tuple{Float64, Float64}=(-10.0, 10.0),\n    θSol::Vector{Float64}=10.0 * randn(Float64, n),\n    std::Float64=200.0, out_times::Float64=7.0)\n\nRandom generate a fitting one-dimensional data problem, storing the data in matrix data and the outliers in vector v.\n\nThis function receives a model(x, θ) function, the number of parameters n, the number of points np to be generated and the number of trusted points p. \n\nIf the n-dimensional vector θSol is provided, then the exact solution will not be random generated. The interval [xMin, xMax] (deprecated) or x_interval for generating the values to evaluate model can also be provided.\n\nIt returns a tuple (data, θSol, outliers) where\n\ndata: (np x 3) array, where each row contains x and model(x, θSol).\nθSol: n-dimensional vector with the exact solution.\noutliers: the outliers of this data set\n\n\n\n\n\n"
},

{
    "location": "api/#RAFF.generate_noisy_data",
    "page": "API",
    "title": "RAFF.generate_noisy_data",
    "category": "function",
    "text": "generate_noisy_data(model::Function, n::Int, np::Int, p::Int;\n    x_interval::Tuple{Float64, Float64}=(-10.0, 10.0),\n    θSol::Vector{Float64}=10.0 * randn(Float64, n),\n    std::Float64=200.0, out_times::Float64=7.0)\n\ngenerate_noisy_data(model::Function, n::Int, np::Int, p::Int,\n    x_interval::Tuple{Float64, Float64})\n\ngenerate_noisy_data(model::Function, n::Int, np::Int, p::Int,\n    θSol::Vector{Float64}, x_interval::Tuple{Float64, Float64})\n\nRandom generate a fitting one-dimensional data problem.\n\nThis function receives a model(x, θ) function, the number of parameters n, the number of points np to be generated and the number of trusted points p. \n\nIf the n-dimensional vector θSol is provided, then the exact solution will not be random generated. The interval [xMin, xMax] (deprecated) or x_interval for generating the values to evaluate model can also be provided.\n\nIt returns a tuple (data, θSol, outliers) where\n\ndata: (np x 3) array, where each row contains x and model(x, θSol).\nθSol: n-dimensional vector with the exact solution.\noutliers: the outliers of this data set\n\n\n\n\n\n"
},

{
    "location": "api/#RAFF.generate_clustered_noisy_data!",
    "page": "API",
    "title": "RAFF.generate_clustered_noisy_data!",
    "category": "function",
    "text": "generate_clustered_noisy_data!(data::Array{Float64, 2},\n    v::Vector{Int}, model::Function, n::Int, np::Int, p::Int,\n    x_interval::Tuple{Float64,Float64},\n    cluster_interval::Tuple{Float64, Float64}; kwargs...)\n\nGenerate a test set with clustered outliers. This version overwrites the content of (np x 3) matrix data and vector v with integer indices to the position of outliers in data.\n\nThe arguments and optional arguments are the same for generate_noisy_data!, with exception of tuple cluster_interval which is the interval to generate the clustered outliers.\n\nIt returns a tuple (data, θSol, outliers) where\n\ndata: (np x 3) array, where each row contains x and model(x, θSol). The same array given as argument\nθSol: n-dimensional vector with the exact solution.\noutliers: the outliers of this data set. The same vector given as argument.\n\n\n\n\n\n"
},

{
    "location": "api/#RAFF.generate_clustered_noisy_data",
    "page": "API",
    "title": "RAFF.generate_clustered_noisy_data",
    "category": "function",
    "text": "generate_clustered_noisy_data(model::Function, n::Int, np::Int,\n    p::Int, x_interval::Tuple{Float64,Float64},\n    cluster_interval::Tuple{Float64, Float64}; kwargs...)\n\ngenerate_clustered_noisy_data(model::Function, n::Int,\n    np::Int, p::Int, θSol::Vector{Float64},\n    x_interval::Tuple{Float64,Float64},\n    cluster_interval::Tuple{Float64, Float64}; kwargs...)\n\nGenerate a test set with clustered outliers.\n\nThe arguments and optional arguments are the same for generate_noisy_data!, with exception of tuple cluster_interval which is the interval to generate the clustered outliers.\n\nIt returns a tuple (data, θSol, outliers) where\n\ndata: (np x 3) array, where each row contains x and model(x, θSol). The same array given as argument\nθSol: n-dimensional vector with the exact solution.\noutliers: the outliers of this data set. The same vector given as argument.\n\n\n\n\n\n"
},

{
    "location": "api/#RAFF.generate_circle",
    "page": "API",
    "title": "RAFF.generate_circle",
    "category": "function",
    "text": "generate_circle(dat_filename::String, np::Int, p::Int;\n    std::Float64=0.1, θSol::Vector{Float64}=1.0*randn(Float64, 3),\n    outTimes::Float64=3.0, interval=(rand(i)*2.0*π for i = 1:np))\n\nGenerate perturbed points in a circle given by θSol and save to dat_filename in RAFF format. Return the np x 4 matrix with data (the 4th column is 0 if the point is \"correct\") and a np - p integer vector containing the points selected to be outliers.\n\ndat_filename is a String with the name of the file to store generated data.\nnp is the number of points to be generated.\np is the number of trusted points to be used in the LOVO approach.\n\nAdditional configuration parameters are\n\nstd: standard deviation.\nθSol: true solution, used for generating perturbed points.\nout_times: deviation for outliers will be out_times * std.\ninterval: any iterable object containing np numbers between 0 and 2π.\n\n\n\n\n\n"
},

{
    "location": "api/#RAFF.generate_ncircle",
    "page": "API",
    "title": "RAFF.generate_ncircle",
    "category": "function",
    "text": "generate_ncircle(dat_filename::String,np::Int, p::Int;\n  std::Float64=0.1, θSol::Vector{Float64}=10.0*randn(Float64, 3),\n  interval=(rand()*2.0*π for i = 1:np))\n\nGenerate perturbed points and uniform noise in a square containing the circle given by θSol and save data to dat_filename in RAFF format. Return the np x 4 matrix with data (the 4th column is 0 if the point is \"correct\") and a np - p integer vector containing the points selected to be outliers.\n\ndat_filename is a String with the name of the file to store generated data.\nnp is the number of points to be generated.\np is the number of trusted points to be used in the LOVO approach.\n\nAdditional configuration parameters are\n\nstd: standard deviation.\nθSol: true solution, used for generating perturbed points.\ninterval: any iterable object containing np numbers between 0 and 2π.\nleftd: number of times the radius of the circle that will be used for computing the lower left corner of the square for generation of the random noise\nlngth: number of times the radius of the circle that will be used for computing the side of the square for generation of the random noise\n\n\n\n\n\n"
},

{
    "location": "api/#RAFF.generate_image_circle",
    "page": "API",
    "title": "RAFF.generate_image_circle",
    "category": "function",
    "text": "generate_image_circle(dat_filename::String, w::Int, h::Int,\n    np::Int, p::Int; std=0.1,\n    θSol::Vector{Float64}=10.0*randn(Float64, 3),\n    interval=(rand()*2.0*π for i = 1:p), thck::Int=2,\n    funcsize=min(w, h))\n\nGenerate perturbed points and uniform noise in a wxh image containing the circle given by θSol and save data to dat_filename in RAFF format. Return the 0-1 matrix representing the black and white image generate.\n\ndat_filename is a String with the name of the file to store generated data.\nw and h are the dimensions of the image\nnp is the number of points to be generated.\np is the number of trusted points to be used in the LOVO approach.\n\nAdditional configuration parameters are\n\nstd: standard deviation.\nθSol: true solution, used for generating perturbed points.\ninterval: any iterable object containing np numbers between 0 and 2π.\nthck: thickness of the point in the image\nfuncsize: size (in pixels) that the function will use in the image.\n\n\n\n\n\n"
},

{
    "location": "api/#RAFF.generate_image_noisy_data",
    "page": "API",
    "title": "RAFF.generate_image_noisy_data",
    "category": "function",
    "text": "function generate_image_noisy_data(dat_filename::String,\nw::Int, h::Int, model::Function, n::Int, np::Int, p::Int;\nx_interval::Tuple{Number, Number}=(-10.0, 10.0),\nθSol::Vector{Float64}=10.0 * randn(Float64, n), std=2,\nthck::Int=2, funcsize=min(w, h))\n\nCreate a file dat_filename with data information to detect model in a wxh image containing random uniform noise. Attention: this function only works with 1-dimensional models.\n\nReturn a black and white matrix representing the image.\n\nThe parameters are\n\ndat_filename: name of the file to save data\nw and h: dimension of the image\nmodel: real-valued model given by a function model(x, θ)\nn: dimension of the parameters of the model\nnp: number of points to be generated\np: number of trusted points that will define the correct points in the model\n\nThe function also accepts the following optional arguments:\n\nx_interval: tuple representing the interval for the x variable\nθSol: vector with the \'exact\' parameters of the solution\nstd: error that will be added to the simulated \'correct\' points\nthck: thickness of the point in the image\nfuncsize: size (in pixels) that the function will use in the image.\n\n\n\n\n\n"
},

{
    "location": "api/#RAFF.model_list",
    "page": "API",
    "title": "RAFF.model_list",
    "category": "constant",
    "text": "This dictionary represents the list of models used in the generation of random tests. Return the tuple (n, model, model_str), where\n\nn is the number of parameters of the model\nmodel is the model of the form m(x, θ), where x are the variables and θ are the parameters\nmodel_str is the string representing the model, used to build random generated problems\n\n\n\n\n\n"
},

{
    "location": "api/#Random-generation-1",
    "page": "API",
    "title": "Random generation",
    "category": "section",
    "text": "RAFF.generate_test_problems\nRAFF.get_unique_random_points\nRAFF.get_unique_random_points!\nRAFF.generate_noisy_data!\nRAFF.generate_noisy_data\nRAFF.generate_clustered_noisy_data!\nRAFF.generate_clustered_noisy_data\nRAFF.generate_circle\nRAFF.generate_ncircle\nRAFF.generate_image_circle\nRAFF.generate_image_noisy_data\nRAFF.model_list"
},

{
    "location": "api/#RAFF.RAFFOutput",
    "page": "API",
    "title": "RAFF.RAFFOutput",
    "category": "type",
    "text": "This type defines the output file for the RAFF algorithm.\n\nRAFFOutput(status::Int, solution::Vector{Float64}, iter::Int,\n           p::Int, f::Float64, nf::Int, nj::Int, outliers::Vector{Int})\n\nwhere\n\nstatus: is 1 if converged and 0 if not\nsolution: vector with the parameters of the model\niter: number of iterations up to convergence\np: number of trusted points\nf: the residual value\nnf: number of function evaluations\nnj: number of Jacobian evaluations\noutliers: the possible outliers detected by the method, for the given p\nRAFFOutput()\n\nCreates a null version of output, equivalent to RAFFOutput(0, [], -1, 0, Inf, -1, -1, [])\n\nRAFFOuput(p::Int)\nRAFFOuput(sol::Vector{Float64}, p::Int)\n\nCreates a null version of output for the given p and a null version with the given solution, respectively.\n\n\n\n\n\n"
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
    "text": "In order to develop a robust code, test problems regarding data fitting with outliers have been created. All the problems are available in directory test/test_problems. All the problems have the formbeginarrayccccc\n	k     \n	x_11  x_12    x_1k  f_1\n	     \n	x_np1  x_np2    x_npk  f_np\nendarraywhere k is the number of variables of the model (not the number of parameters!), np is the number of data points to be adjusted, x_ij are the values selected for parameter j in experiment i and f_i is the result obtained for experiment i."
},

{
    "location": "advanced/#Script-files-1",
    "page": "Advanced",
    "title": "Script files",
    "category": "section",
    "text": "During the development and testing of RAFF several scripts and pieces of Julia code have been created. Those files are mostly related to the automated generation of test problems and visualization of the solutions. All those files are located in the test/scripts directory.We explain each file below, so maybe more advanced users can modify and re-use the code to generate their own problems.calc_ratio.jl: this script contains a function that generates several tests for the same model and solve each one with RAFF. Then it prints the ratio of outliers that have successfully been detected but the method.\ndraw.jl: This script draws the solution of a given problem and also its data, which was taken from a file. Examples of such files are located in test/test_problems.\nrun_raff.jl: this script simply loads a problem data from a file, selects a given model and runs RAFF. Examples of such files are located in test/test_problems.\nrun_lmlovo.jl: the same as before, but only for lmlovo function. Used mostly for testing.\ngenerate_fit_tests.jl: script for generating random test problem files, using the pre-defined models given by RAFF.model_list. This function cannot be called inside Julia, since it uses ArgParse package.\ngen_circle.jl: specific script for generating random test problems related to the detection of circles in the plane. It also provides functions to draw the problem and the solution, which differ from the draw.jl script above.\nrun_performance_tests.jl: script for generating some performance tests, so we can compare different versions of RAFF."
},

{
    "location": "randomgeneration/#",
    "page": "Random generation",
    "title": "Random generation",
    "category": "page",
    "text": ""
},

{
    "location": "randomgeneration/#Random-generation-of-problems-1",
    "page": "Random generation",
    "title": "Random generation of problems",
    "category": "section",
    "text": "RAFF.jl contains several methods for the generation of artificial datasets, in order to simulate noisy data with outliers. See the API for details about the possibilities of random generation of data with outliers."
},

{
    "location": "randomgeneration/#Simple-example-the-exponential-function-1",
    "page": "Random generation",
    "title": "Simple example - the exponential function",
    "category": "section",
    "text": "First, it is necessary to load RAFF.jl and the desired model to generate the problem.using RAFF\n\nn, model, = RAFF.model_list[\"expon\"]If an exact solution is not provided, RAFF.jl will generate a random one, but this can destroy the shape of some models. Therefore, in this example, we will provide a hint for a nice exponential model. In this example, we will generate 20 random points with 2 outliers in the interval 1 30.exact_sol = [5000.0, 4000.0, 0.2]\n\ninterv = (1.0, 30.0)\n\nnp = 20\n\np = 18Before calling the generating function, we fix the random seed, so this example is the same everywhere.using Random\n\nRandom.seed!(12345678)\n\ndata, = generate_noisy_data(model, n, np, p, x_interval=interv, θSol=exact_sol)Now we use the script test/script/draw.jl (see Advanced section) to visualize the data generated. draw.jl uses the PyPlot.jl package, so it is necessary to install it before running this example.julia> include(\"test/scripts/draw.jl\")\n\njulia> draw_problem(data, model_str=\"expon\")If you are interested in running RAFF, just call the raff method. Attention: we have to drop the last column of data, since it contains information regarding the noise added to the outliers.r = raff(model, data[:, 1:end - 1], n, MAXMS=10)If all the steps were successful, after commandjulia> draw_problem(data, model_str=\"expon\", raff_output=r)the following picture should appear(Image: Exponential example)"
},

]}
