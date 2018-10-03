# This file contains utility functions to deal with distributed
# computing

"""

    update_best(channel::RemoteChannel, bestx::SharedArray{Float64, 1})

Listen to a `channel` for results found by LMlovo. If there is an
improvement for the objective function, the shared array `bestx` is
updated.

**Atention**: There might be an unstable state if there is a process
  reading `bestx` while this function is updating it. This should not
  be a problem, since it is used as a starting point.

"""
function update_best(channel::RemoteChannel, bestx::SharedArray{Float64, 1})

    @debug("Running updater.")

    n = length(bestx)

    N::Int = 0
    
    while isopen(channel)

        x = try
            
            take!(channel)
            
        catch e

            if isa(e, InvalidStateException)
                break
            end
            
            @warn("Something wrong when reading from channel. Will skip.", e)
            
            continue
            
        end
        
        @debug("Updater has read values from channel.")

        for i = 1:n
        
            bestx[i] = (N * bestx[i] + x[i]) / (N + 1)

        end

        N += 1
        
    end

    @debug("Update channel closed. Exiting thread.")
    
end

"""

    consume_tqueue(bqueue::RemoteChannel, tqueue::RemoteChannel,
                   bestx::SharedArray{Float64, 1},
                   v::SharedArray{Float64, 2},
                   vs::SharedArray{Int,1},
                   vf::SharedArray{Float64, 1}, model::Function,
                   gmodel!::Function, data::Array{Float64, 2}, n::Int,
                   pliminf::Int, plimsup::Int, MAXMS::Int,
                   seedMS::MersenneTwister)

This function represents one worker, which runs LMlovo in a multistart
fashion.

It takes a job from the RemoteChannel `tqueue` and runs `LMlovo`
function to it. Saves the best results found to the shared arrays `v`
(best solution), `vs` (convergence status) and `vf` (objective
function value at the best solution). All the other arguments are the
same for `praff` function.

"""
function consume_tqueue(bqueue::RemoteChannel, tqueue::RemoteChannel,
                        bestx::SharedArray{Float64, 1},
                        v::SharedArray{Float64, 2},
                        vs::SharedArray{Int, 1},
                        vf::SharedArray{Float64, 1}, model::Function,
                        gmodel!::Function, data::Array{Float64, 2},
                        n::Int, pliminf::Int, plimsup::Int,
                        MAXMS::Int, seedMS::MersenneTwister)

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

        if (p < pliminf) || (p > plimsup)

            @warn("Invalid value for task: $(p). Will skip task.")

            continue

        end
        
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
                if s == 1

                    try
                        
                        put!(bqueue, x)
                        
                        @debug("Added new point to queue.", x, f)

                    catch e

                        @warn(string("Problems when saving best point found in queue. ",
                                     "Will skip this step"), e)

                    end

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

        @debug("Finished. p = $(p) and f = $(wbestf). -> $(wbestx)")

    end

end
