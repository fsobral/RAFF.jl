using Documenter, DocumenterTools,RAFF

makedocs(
	assets = ["assets/favicon.ico"],
	sitename = "RAFF- Robust Algebraic Fitting Function",
	pages = ["Overview" => "index.md",
	"Getting Started" => "gettingstarted.md",
	"API" => "summary.md"],
	#html_prettyurls = false
	format = Documenter.HTML(prettyurls = false)
	)
	
