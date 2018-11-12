using Documenter, RAFF

makedocs(format = :html,
	 assets = ["assets/favicon.ico"],
	 sitename = "RAFF- Robust Algebraic Fitting Function",
	 pages = ["Overview" => "index.md",
		  "Getting Started" => "gettingstarted.md"]
	)
	
