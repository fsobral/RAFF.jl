# This script generates and draws examples for the circle detection

using PyPlot
using Printf
using DelimitedFiles
using RAFF
using FileIO, Images, ImageView

"""

    gen_circle(np::Int, p::Int; std::Float64=0.1,
               θSol::Vector{Float64}=10.0*randn(Float64, 3),
               outTimes::Float64=5.0,
               interval::Vector{Float64}=rand(np)*2.0*π)

Generate perturbed points in a circle given by `θSol`. Construct a
test file for `RAFF`.

"""
function gen_circle(np::Int, p::Int; std::Float64=0.1,
                    θSol::Vector{Float64}=1.0*randn(Float64, 3),
                    outTimes::Float64=3.0, interval::Vector{Float64}=rand(np)*2.0*π)

    ρ = (α, ρ) -> [ρ * cos(α) + θSol[1], ρ * sin(α) + θSol[2]]
    f = (x) -> (x[1] - θSol[1])^2 + (x[2] - θSol[2])^2 - θSol[3]^2

    data = Array{Float64, 2}(undef, np, 4)

    # Points selected to be outliers
    v = RAFF.get_unique_random_points(np, np - p)

    for (i, α) in enumerate(interval)
    
        pt = ρ(α, θSol[3] + std * randn())

        data[i, 3:4] .= 0.0 #f(pt)

        # Follow the noise idea of J. Yu, H. Zheng, S. R. Kulkarni,
        # and H. V. Poor, "Two-Stage Outlier Elimination for Robust
        # Curve and Surface Fitting," EURASIP J. Adv. Signal Process.,
        # vol. 2010, no. 1, p. 154891, Dec. 2010.
        if i in v

            # pt = ρ(α, θSol[3] * (1.0 + 2 * rand()) * outTimes * std * sign(randn()))

            pt[1] += outTimes * std * randn()
            pt[2] += outTimes * std * randn()
            
            data[i, 4] = 1.0
            
        end
        
        data[i, 1:2] = pt

    end

    open("/tmp/output.txt", "w") do fp
    
        # Dimension of the domain of the function to fit
        @printf(fp, "%d\n", 2)

        for k = 1:np

            @printf(fp, "%20.15f %20.15f %20.15f %1d\n",
                    data[k, 1], data[k, 2], data[k, 3], Int(k in v))

        end

    end

    return data, v

end

function gen_ncircle(np::Int, p::Int; std::Float64=0.1,
                    θSol::Vector{Float64}=10.0*randn(Float64, 3),
                    interval::Vector{Float64}=rand(p)*2.0*π)

    ρ = (α, ρ) -> [ρ * cos(α) + θSol[1], ρ * sin(α) + θSol[2]]
    f = (x) -> (x[1] - θSol[1])^2 + (x[2] - θSol[2])^2 - θSol[3]^2

    data = Array{Float64, 2}(undef, np, 4)

    for (i, α) in enumerate(interval)
    
        pt = ρ(α, θSol[3] + std * randn())

        data[i, 1:2]  = pt
        data[i, 3:4] .= 0.0 #f(pt)

    end

    # Just random noise

    v = Vector{Int}(undef, np - p)
    
    for i = p + 1:np

        data[i, 1] = θSol[1] - 2.0 * θSol[3] + rand() * 4.0 * θSol[3]
        data[i, 2] = θSol[2] - 2.0 * θSol[3] + rand() * 4.0 * θSol[3]
        data[i, 3] = 0.0
        data[i, 4] = 1.0

        v[i - p] = i

    end

    open("/tmp/output.txt", "w") do fp
    
        # Dimension of the domain of the function to fit
        @printf(fp, "%d\n", 2)

        for k = 1:np

            @printf(fp, "%20.15f %20.15f %20.15f %1d\n",
                    data[k, 1], data[k, 2], data[k, 3], Int(k in v))

        end

    end

    return data, v

end


"""

Draw the points generated by the previous function.

"""
function draw_circle(data, outliers)

    np, = size(data)
    
    c = zeros(np)
    
    c[outliers] .= 1.0
    
    PyPlot.scatter(data[:, 1], data[:, 2], c=c, marker="o", s=50.0, linewidths=0.2,
                   cmap=PyPlot.cm["Paired"], alpha=0.9)
    
    PyPlot.axis("scaled")

    PyPlot.xticks([])

    PyPlot.yticks([])

    PyPlot.savefig("/tmp/circle.png", dpi=150, bbox_inches="tight")

end

"""

Draw the points and the solutions obtained. Save the picture in a file.

"""
function draw_circle_sol(M; model_str="circle", raff_output=nothing, other_sols...)

    x  = M[:, 1]
    y  = M[:, 2]
    co = M[:, 4]

    t = 0:0.1:2.1 * π
    
    ptx = (α, ρ, d) -> ρ * cos(α) + d[1]
    pty = (α, ρ, d) -> ρ * sin(α) + d[2]    

    # Plot data

    true_outliers = findall(co .!= 0.0)

    PyPlot.scatter(x[co .== 0.0], y[co .== 0.0], color=PyPlot.cm."Pastel1"(2.0/9.0),
                   marker="o", s=50.0, linewidths=0.2)

    PyPlot.scatter(x[co .!= 0.0], y[co .!= 0.0], color=PyPlot.cm."Pastel1"(2.0/9.0),
                   marker="^", s=25.0, linewidths=0.2, label="Outliers")

    if raff_output != nothing
        
        n, model, modelstr = RAFF.model_list[model_str]

        fSol = raff_output.solution
        
        modl1x = (α) -> ptx(α, fSol[3], fSol[1:2])
        modl1y = (α) -> pty(α, fSol[3], fSol[1:2])

        PyPlot.plot(modl1x.(t), modl1y.(t), color=PyPlot.cm."Set1"(2.0/9.0))

        # Draw outliers found by RAFF
        
        true_positives = intersect(true_outliers, raff_output.outliers)
        false_positives = setdiff(raff_output.outliers, true_positives)

        if length(false_positives) > 0
        
            PyPlot.scatter(x[false_positives], y[false_positives],
                           color=PyPlot.cm."Pastel1"(0.0/9.0), marker="o",
                           linewidths=0.2, edgecolors="k", s=50.0, label="False positives")

        end

        if length(true_positives) > 0
        
            PyPlot.scatter(x[true_positives], y[true_positives],
                           color=PyPlot.cm."Pastel1"(0.0/9.0), marker="^",
                           s=50.0, linewidths=0.2, edgecolors="k", label="Identified outliers")

        end

    end

    PyPlot.legend(loc="best")

    PyPlot.axis("scaled")

    PyPlot.xticks([])

    PyPlot.yticks([])

    PyPlot.savefig("/tmp/circle.png", dpi=150, bbox_inches="tight")
    
end

function draw_circle_sol(tSol, fSol, lsSol)

    datafile = "/tmp/output.txt"

    fp = open(datafile, "r")

    N = parse(Int, readline(fp))

    M = readdlm(fp)

    close(fp)

    x  = M[:, 1]
    y  = M[:, 2]
    ρ  = M[:, 3]
    co = M[:, 4]

    t = [0:0.1:2.1 * π;]
    
    ptx = (α, ρ, d) -> ρ * cos(α) + d[1]
    pty = (α, ρ, d) -> ρ * sin(α) + d[2]

    # True solution
    
    pptx = (α) -> ptx(α, tSol[3], tSol[1:2])
    ppty = (α) -> pty(α, tSol[3], tSol[1:2])

    PyPlot.plot(pptx.(t), ppty.(t), "b--", label="True solution")
    
    # RAFF solution
    
    pptx = (α) -> ptx(α, fSol[3], fSol[1:2])
    ppty = (α) -> pty(α, fSol[3], fSol[1:2])

    PyPlot.plot(pptx.(t), ppty.(t), "g-", label="RAFF")
    
    # # LS solution
    
    # pptx = (α) -> ptx(α, lsSol[3], lsSol[1:2])
    # ppty = (α) -> pty(α, lsSol[3], lsSol[1:2])

    # PyPlot.plot(pptx.(t), ppty.(t), "r-", label="Least squares")

    PyPlot.scatter(x[co .== 0.0], y[co .== 0.0], color=PyPlot.cm."Pastel1"(2.0/9.0),
                   marker="o", s=50.0, linewidths=0.2)

    PyPlot.scatter(x[co .!= 0.0], y[co .!= 0.0], color=PyPlot.cm."Pastel1"(1.0/9.0),
                   marker=".", s=25.0, linewidths=0.2, label="Outliers")

    PyPlot.legend(loc=4)
    
    PyPlot.axis("scaled")

    PyPlot.xticks([])

    PyPlot.yticks([])

    PyPlot.savefig("/tmp/circle.png", dpi=150, bbox_inches="tight")

end

"""

    draw_img_sol(imgfile, sol; thck::Int=2)

Draw the image and the solution found, using JuliaImage
packages. `imagefile` is the path to the image and `sol` is a
3-dimensional vector with the solution.

Extra parameter `thck` defines the thickness of the circle.

"""
function draw_img_sol(imgfile, sol; thck::Int=2)

    img = load(imgfile)

    h, w = size(img)
    
    t = [0:0.01:2.1 * π;]
    
    ptx = (α, ρ, d) -> ρ * cos(α) + d[1]
    pty = (α, ρ, d) -> ρ * sin(α) + d[2]

    # Solution
    
    pptx = (α) -> ptx(α, sol[3], sol[1:2])
    ppty = (α) -> pty(α, sol[3], sol[1:2])

    for α in t

        i = Int(round(pptx(α)))
        j = Int(round(ppty(α)))

        !((1 <= i <= h) && (1 <= j <= w)) && continue
        
        img[max(1, i - thck):min(i + thck, h), max(1, j - thck):min(j + thck, w)] .= RGB(1.0, 0.0, 0.0)

    end
    
    ImageView.imshow(img)
    
end
