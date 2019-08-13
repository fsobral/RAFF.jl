using PyPlot,LinearAlgebra,Random
"""
```
gen_ellipse(F1:Vector{Float64},F2=Vector{Float64},
	a::Float64,p::Int,npun::Int,ntot::Int,
	tmin::Float64,tmax::Float64)
```

Generate perturbed points in a ellipse given focus F1 and F2 major axis a. It is mandatory to define the number of points `ntot` the of exacts points `p` and the number of point to compose the ellipse `npun`. 

`ntot-npun` is the number of random points
`npun-p` is the number of points of the ellipse which are normally-distributed.
`tmin` and `tmax` are parameter to define limits for the figure. 

## Example

```julia-repl
julia> gen_ellipse([1.0,2.0],[3.0,5.0],5.0,0,100,1000,-10,10)
```

"""
function gen_ellipse(F1,F2,a,p,npun,ntot,tmin,tmax)
	c = norm(F2-F1,2)
	if a<=c 
		error("a<=c, it isn't an elipse")
	end
	b = sqrt(a^2-c^2)
	# let us assume β is the rotation angle from (1,0)
	cosβ = (F2-F1)[1]/norm(F2-F1,2)
	sinβ = sqrt(1.0-cosβ)
	
	center = 0.5*(F1+F2)
	rng = MersenneTwister(1234)
	x = zeros(ntot)
	y = zeros(ntot)
	θ = shuffle(rng,[0:(2*π)/npun:2*π;])
	
	for i=1:p
		x[i] = center[1]*cosβ+b*cos(θ[i])*cosβ+center[2]*sinβ+a*sin(θ[i])*sinβ
		y[i] = -center[1]*sinβ -b*cos(θ[i])*sinβ+center[2]*cosβ+a*sin(θ[i])*cosβ
	end
	for i=p+1:npun
		anoise = a+0.1*randn()
		bnoise = b+0.1*randn()
		x[i] = center[1]*cosβ+bnoise*cos(θ[i])*cosβ+center[2]*sinβ+anoise*sin(θ[i])*sinβ
		y[i] = -center[1]*sinβ -bnoise*cos(θ[i])*sinβ+center[2]*cosβ+anoise*sin(θ[i])*cosβ
	end

	x[npun+1:end]= tmin.+(tmax-tmin)*rand(ntot-npun)
	y[npun+1:end]= tmin.+(tmax-tmin)*rand(ntot-npun)
	plot(x,y,".")

	return [x y]
end

