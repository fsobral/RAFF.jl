# This example shows how to get mouse clicks in Julia

# Packages needed
#
# add Gtk
# add GtkReactive
# add Graphics
# add Colors

using Gtk, Gtk.ShortNames, GtkReactive, Graphics, Colors

win = Window("Drawing")
c = canvas(UserUnit)       # create a canvas with user-specified coordinates
push!(win, c)

const lines = Signal([])   # the list of lines that we'll draw
const newline = Signal([]) # the in-progress line (will be added to list above)

const drawing = Signal(false)  # this will become true if we're actively dragging

# c.mouse.buttonpress is a `Reactive.Signal` that updates whenever the
# user clicks the mouse inside the canvas. The value of this signal is
# a MouseButton which contains position and other information.

# We're going to define a callback function that runs whenever the
# button is clicked. If we just wanted to print the value of the
# returned button object, we could just say
#     map(println, c.mouse.buttonpress)
# However, here our function is longer than `println`, so
# we're going to use Julia's do-block syntax to define the function:
sigstart = map(c.mouse.buttonpress) do btn
    # This is the beginning of the function body, operating on the argument `btn`
    if btn.button == 1 && btn.modifiers == 0 # is it the left button, and no shift/ctrl/alt keys pressed?
        push!(drawing, true)   # activate dragging
        push!(newline, [btn.position])  # initialize the line with the current position
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
        # append our new line to the overall list
        push!(lines, push!(value(lines), value(newline)))
        # For the next click, make sure `newline` starts out empty
        push!(newline, [])
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

# Because `draw` isn't a one-line function, we again use do-block syntax:
redraw = draw(c, lines, newline) do cnvs, lns, newl  # the function body takes 3 arguments
    fill!(cnvs, colorant"white")   # set the background to white
    set_coordinates(cnvs, BoundingBox(0, 1, 0, 1))  # set coordinates to 0..1 along each axis
    ctx = getgc(cnvs)   # gets the "graphics context" object (see Cairo/Gtk)
    for l in lns
        drawline(ctx, l, colorant"blue")  # draw old lines in blue
    end
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

