using Documenter, DocumenterTools,RAFF

makedocs(
	assets = ["assets/favicon.ico"],
	sitename = "RAFF- Robust Algebraic Fitting Function",
	pages = ["Overview" => "index.md",
	"Tutorial"=> "tutorial.md",
	"API" => "api.md"],
	#html_prettyurls = false
	format = Documenter.HTML(prettyurls = false)
	)
	
