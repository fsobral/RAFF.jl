# This file contains utility functions to deal with distributed
# computing

"""

    update_best(channel::RemoteChannel, bestx::SharedArray{Float64, 1})

Listen to a `channel` for results found by lmlovo. If there is an
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

This function represents one worker, which runs lmlovo in a multistart
fashion.

It takes a job from the RemoteChannel `tqueue` and runs `lmlovo`
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

    # my_bestx::Vector{Float64} = zeros(Float64, n)
    
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

        if (p.start < pliminf) || (p.stop > plimsup) ||
           (length(p) == 0)

            @warn("Invalid value for task: $(p). Will skip task.")

            continue

        end

        for k in p
        
            # Starting point
            wbestx = zeros(Float64, n)
            wbestf = Inf
            ws = 0
            
            # Multi-start strategy
            for j = 1:MAXMS

                # New random starting point
                # x = randn(seedMS, n)
                # x .= 10.0 .* x .+ my_bestx
                # Best results with this one
                x = zeros(Float64, n)
                
                # Call function and store results
                s, x, iter, p_, f = lmlovo(model, gmodel!, x, data, n, k)
                
                if f < wbestf
                    
                    # Send asynchronously the result to channel if success
                    if s == 1

                        @async try
                            
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

            # @async my_bestx .= bestx

            @async begin
                
                ind = k - pliminf + 1
                
                v[:, ind] .= wbestx
                vf[ind]    = wbestf
                vs[ind]    = ws

            end

            @debug("Finished. p = $(k) and f = $(wbestf). -> $(wbestx)")

        end

    end

end

"""

    check_and_close(bqueue::RemoteChannel, tqueue::RemoteChannel,
                    futures::Vector{Future}; secs::Float64=0.1)

Check if there is at least one worker process in the vector of
`futures` that has not prematurely finished. If there is no alive
worker, close task and best queues `tqueue` and `bqueue`,
respectively.

"""
function check_and_close(bqueue::RemoteChannel, tqueue::RemoteChannel,
                         futures::Vector{Future}; secs::Float64=0.1)

    n_alive = length(futures)

    @debug("Checking worker status.")

    for (i, f) in enumerate(futures)

        if timedwait(()->isready(f), secs) == :ok

            @warn("Worker $(i) has finished prematurely.")

            n_alive -= 1

        end

    end

    @debug("Workers online: $(n_alive)")
    
    if n_alive == 0

        @warn("No live worker found. Will close queues and finish.")

        close(bqueue)

        close(tqueue)

    end

end
