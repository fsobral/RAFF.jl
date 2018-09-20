# This example shows a simple use for a channel and a feeder

using Distributed
using Logging
using Base.CoreLogging

# Set Debug for Logging
global_logger(ConsoleLogger(stdout, Logging.Debug))

function read_and_print(channel::Channel{Tuple{Vector{Int}, Int}})

    @debug("Running updater.")
    
    while true

        x, i = try
            take!(channel);
        catch e

            if isa(e, InvalidStateException)
                break
            end
            
            @warn("Something wrong when reading from channel. Will skip.", e)
            
            continue
            
        end

        @debug("Updater has read values from channel. x = $(x) and i = $(i).")
       
    end

    @debug("Update channel closed. Exiting thread.")
    
end

c = Channel{Tuple{Vector{Int}, Int}}(0)

# Task waits for the channel to be fed
task1 = @async read_and_print(c)

# Feeding the channel
task2 = @async for i = 1:10
    put!(c, (ones(Int, 3), i))
end

sleep(1)

# Wrong feeding
put!(c, 1)

# This should close `task1`
close(c)

# Just to be sure
sleep(1)

println("Task1 done = $(istaskdone(task1))")

println("Task2 done = $(istaskdone(task2))")
