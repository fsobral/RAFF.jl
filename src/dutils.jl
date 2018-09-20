# This file contains utility functions to deal with distributed
# computing

"""

    update_best()

Listen to a `channel` for results found by LMlovo. If there is an
improvement for the objective function, the shared array `bestx` is
updated.

**Atention**: There might be an unstable state if there is a process
  reading `bestx` while this function is updating it.

"""
function update_best(channel::RemoteChannel, bestx::SharedArray{Float64, 1})

    @debug("Running updater.")

    n = length(bestx)

    N::Int = 0
    
    while isopen(channel)

        (x, f) = try
            
            take!(channel)
            
        catch e

            if isa(e, InvalidStateException)
                break
            end
            
            @warn("Something wrong when reading from channel. Will skip.", e)
            
            continue
            
        end
        
        @debug("Updater has read values from channel. f = $(f).")

        for i = 1:n
        
            bestx[i] = (N * bestx[i] + x[i]) / (N + 1)

        end

        N += 1
        
    end

    @debug("Update channel closed. Exiting thread.")
    
end
