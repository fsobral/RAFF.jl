using Documenter, DocumenterTools,RAFF

makedocs(
	assets = ["assets/favicon.ico"],
	sitename = "RAFF- Robust Algebraic Fitting Function",
	pages = ["Overview" => "index.md",
	"Getting Started" => "gettingstarted.md",
	"API" => "api.md"],
	#html_prettyurls = false
	format = Documenter.HTML(prettyurls = false)
	)
	
