# Examples

In addition to the examples given in the [Tutorial](@ref), the
`examples/` directory contains another ways of using RAFF. Currently,
we only provide an example on how to load a problem from file, solve
it using `RAFF` and visually check the results.

  - `cubic.jl`: this example solves a problem using a cubic model,
    with 4 parameters. The example also illustrates how to use
    [`RAFF.model_list`](@ref) utility structure in order to load
    pre-defined models.
	
  - `draw_and_detect.jl`: this nice example uses
    [`GtkReactive.jl`](https://github.com/JuliaGizmos/GtkReactive.jl)
    to show a graphic application of `RAFF.jl` to the detection of
    circles drawn by the user. The user can also see the difference
    between the LOVO approach and the traditional least squares
    technique.
