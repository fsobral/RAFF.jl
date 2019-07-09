# This example shows how to use RAFF.jl to detect circles drawn from
# the user. This example was inspired in `GtkReactive.jl` drawing
# example.

# Packages needed
#
# add Gtk
# add GtkReactive
# add Graphics
# add Colors

using Gtk, Gtk.ShortNames, GtkReactive, Graphics, Colors
using RAFF

win = Window("Drawing") |> (bx = Box(:v))

set_gtk_property!(bx, :spacing, 10)

cb = GtkReactive.dropdown(["RAFF", "Least Squares"])
push!(bx, cb)

c = canvas(UserUnit, 200, 200)       # create a canvas with user-specified coordinates
push!(bx, c)

const moddraw = Signal([]) # the model points
const newline = Signal([]) # the in-progress line (will be added to list above)

const drawing = Signal(false)  # this will become true if we're actively dragging

choice = map(x->(push!(moddraw, []);value(x)), cb)

function run_raff()

    n, model, = RAFF.model_list["circle"]

    np, = size(value(newline))

    data = Array{Float64, 2}(undef, np, 3)

    l = value(newline)
    
    for i in 1:np

        data[i,1] = l[i].x
        data[i,2] = l[i].y
        data[i,3] = 0.0

    end

    r = if value(choice) == "RAFF"
        
        raff(model, data, n, ftrusted=0.7, MAXMS=10)

    else

        raff(model, data, n, ftrusted=1.0, MAXMS=10)

    end
    println(r)
    # Build points for drawing the circle
    p = (α) -> [abs(r.solution[3]) * cos(α) + r.solution[1], abs(r.solution[3]) * sin(α) + r.solution[2]]
    
    push!(moddraw, [p(i) for i in LinRange(0.0, 2 * π, 100)])
    
end

# If has changed the selection box, then run RAFF again.
sigc = map(choice) do c

    run_raff()
    
end

sigstart = map(c.mouse.buttonpress) do btn
    # This is the beginning of the function body, operating on the argument `btn`
    if btn.button == 1 && btn.modifiers == 0 # is it the left button, and no shift/ctrl/alt keys pressed?
        push!(drawing, true)   # activate dragging
        push!(newline, [btn.position])  # initialize the line with the current position
        push!(moddraw, [])
    end
end

const dummybutton = MouseButton{UserUnit}()
# See the Reactive.jl documentation for `filterwhen`
sigextend = map(filterwhen(drawing, dummybutton, c.mouse.motion)) do btn
    # while dragging, extend `newline` with the most recent point
    push!(newline, push!(value(newline), btn.position))
end

sigend = map(c.mouse.buttonrelease) do btn
    if btn.button == 1
        push!(drawing, false)  # deactivate dragging
        run_raff()
    end
end

function drawline(ctx, l, color)
    isempty(l) && return
    p = first(l)
    move_to(ctx, p.x, p.y)
    set_source(ctx, color)
    for i = 2:length(l)
        p = l[i]
        line_to(ctx, p.x, p.y)
    end
    stroke(ctx)
end

function drawcircle(ctx, l, color)
    isempty(l) && return
    p = first(l)
    move_to(ctx, p[1], p[2])
    set_source(ctx, color)
    for i = 2:length(l)
        p = l[i]
        line_to(ctx, p[1], p[2])
    end
    stroke(ctx)
end
    
# Because `draw` isn't a one-line function, we again use do-block syntax:
redraw = draw(c, moddraw, newline) do cnvs, circ, newl  # the function body takes 3 arguments
    fill!(cnvs, colorant"white")   # set the background to white
    set_coordinates(cnvs, BoundingBox(0, 1, 0, 1))  # set coordinates to 0..1 along each axis
    ctx = getgc(cnvs)   # gets the "graphics context" object (see Cairo/Gtk)
    drawcircle(ctx, circ, colorant"blue")  # draw old lines in blue
    drawline(ctx, newl, colorant"red")    # draw new line in red
end

showall(win)

#If we are not in a REPL
if (!isinteractive())

    # Create a condition object
    c = Condition()

    # Get the window
    # win = guidict["gui"]["window"]

    # Notify the condition object when the window closes
    signal_connect(win, :destroy) do widget
        notify(c)
    end

    # Wait for the notification before proceeding ...
    wait(c)
end

