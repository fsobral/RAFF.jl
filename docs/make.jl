using Documenter, RAFF

push!(LOAD_PATH, "../test/scripts/")

makedocs(
    assets = ["assets/favicon.ico"],
    sitename = "RAFF- Robust Algebraic Fitting Function",
    pages = ["Overview" => "index.md",
	     "Tutorial"=> "tutorial.md",
	     "Examples" => "examples.md",
             "API" => "api.md",
             "Advanced" => "advanced.md",
             "Random generation" => "randomgeneration.md"],
    #html_prettyurls = false
    #format = Documenter.HTML(prettyurls = false),
    modules = [RAFF]
	)
	
deploydocs(
	repo = "github.com/fsobral/RAFF.jl.git"
	)
