using Documenter, RAFF

makedocs(
    sitename = "RAFF- Robust Algebraic Fitting Function",
    pages = ["Overview" => "index.md",
	     "Tutorial"=> "tutorial.md",
	     "Examples" => "examples.md",
             "API" => "api.md",
             "Advanced" => "advanced.md"],
    #html_prettyurls = false
    #format = Documenter.HTML(prettyurls = false),
    format = Documenter.HTML(assets = ["assets/favicon.ico"]),
    modules = [RAFF]
	)
	
deploydocs(
	repo = "github.com/fsobral/RAFF.jl.git"
	)
