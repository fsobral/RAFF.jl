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

"""

This function represents one worker, which runs LMlovo in a multistart
fashion.

"""
function consume_tqueue(bqueue::RemoteChannel, tqueue::RemoteChannel,
                        bestx::SharedArray{Float64, 1},
                        v::SharedArray{Float64, 2},
                        vs::SharedArray{Int, 1},
                        vf::SharedArray{Float64, 1}, model::Function,
                        gmodel!::Function, data::Array{Float64, 2},
                        n::Int, pliminf::Int, MAXMS::Int, seedMS::MersenneTwister)

    @debug("Started worker $(myid())")

    while isopen(tqueue)

        p = try
            
            take!(tqueue)
            
        catch e

            if isa(e, InvalidStateException)
                break
            end
            
            @warn("Something wrong when reading task. Will skip task.", e)
            
            continue
            
        end

        @debug("Received task $(p)")
        
        # Starting point
        wbestx = zeros(Float64, n)
        wbestf = Inf
        ws = 0
    
        # Multi-start strategy
        for j = 1:MAXMS

            # New random starting point
            x = randn(seedMS, n)
            x .= 10.0 .* x .+ bestx
            
            # Call function and store results
            s, x, iter, p_, f = LMlovo(model, gmodel!, x, data, n, p)
            
            if f < wbestf
            
                # Send result to channel if success
                if s == 0

                    put!(bqueue, (x, f))

                    @debug("Added new point to queue.", x, f)

                end
                
                # Save best
                wbestx .= x
                wbestf  = f
                ws      = s
                
            end

        end

        ind = p - pliminf + 1
        
        v[:, ind] .= wbestx
        vf[ind]    = wbestf
        vs[ind]    = ws

        @debug("Finished. p = $(p) and f = $(wbestf).-> $(wbestx)")

    end

end
