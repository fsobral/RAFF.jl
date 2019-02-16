# This file contains utility functions to deal with distributed
# computing

"""
    update_best(channel::RemoteChannel, bestx::SharedArray{Float64, 1})

Listen to a `channel` for results found by lmlovo. If there is an
improvement for the objective function, the shared array `bestx` is
updated.

**Attention**: There might be an unstable state if there is a process
  reading `bestx` while this function is updating it. This should not
  be a problem, since it is used as a starting point.

**Attention 2**: this function is currently out of use.

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

    function consume_tqueue(bqueue::RemoteChannel, tqueue::RemoteChannel,
                            squeue::RemoteChannel, model::Function, gmodel!::Function,
                            data::Array{Float64, 2}, n::Int, pliminf::Int,
                            plimsup::Int, MAXMS::Int, seedMS::MersenneTwister)

This function represents one worker, which runs lmlovo in a multistart
fashion.

It takes a job from the RemoteChannel `tqueue` and runs `lmlovo`
function to it. It might run using a multistart strategy, if
`MAXMS>1`. It sends the best results found for each value obtained in
`tqueue` to channel `squeue`, which will be consumed by the main
process. All the other arguments are the same for [`praff`](@ref)
function.

"""
function consume_tqueue(bqueue::RemoteChannel, tqueue::RemoteChannel,
                        squeue::RemoteChannel,
                        model::Function, gmodel!::Function,
                        data::Array{Float64, 2}, n::Int, pliminf::Int,
                        plimsup::Int, MAXMS::Int,
                        seedMS::MersenneTwister)

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

        if (p.start < pliminf) || (p.stop > plimsup) ||
           (length(p) == 0)

            @warn("Invalid value for task: $(p). Will skip task.")

            continue

        end

        for k in p

            wbest = RAFFOutput(0, [], -1, k, Inf, [])
        
            # Multi-start strategy
            for j = 1:MAXMS

                # New random starting point
                x = randn(seedMS, n)
                
                # Call function and store results
                rout = lmlovo(model, gmodel!, x, data, n, k)

                (rout.status == 1) && (rout.f < wbest.f) && (wbest = rout)

                # This block is related to a strategy of smart
                # starting points for the multistart
                # process. Currently, it makes no sense to use it.
                
                # if rout.f < wbest.f
                    
                #     # Send asynchronously the result to channel if success
                #     if rout.status == 1

                #         @async try
                            
                #             put!(bqueue, rout.solution)
                            
                #             @debug("Added new point to queue.", rout.solution, rout.f)

                #         catch e

                #             @warn(string("Problems when saving best point found in queue. ",
                #                          "Will skip this step"), e)

                #         end

                #     end
                    
                # end

            end

            @debug("Finished. p = $(k) and f = $(wbest.f).")

            try

                put!(squeue, wbest)

            catch e

                if isa(e, InvalidStateException)

                    @warn("Solution queue prematurely closed. Unable to save solution for p=$(k).")

                    return

                end

                @warn("Something wrong when sending the solution to queue for p=$(k).", e)

            end

        end

    end

end

"""

    check_and_close(bqueue::RemoteChannel, tqueue::RemoteChannel,
                    squeue::RemoteChannel, futures::Vector{Future};
                    secs::Float64=0.1)

Check if there is at least one worker process in the vector of
`futures` that has not prematurely finished. If there is no alive
worker, close task, solution and best queues, `tqueue`, `squeue` and
`bqueue`, respectively.

"""
function check_and_close(bqueue::RemoteChannel, tqueue::RemoteChannel,
                         squeue::RemoteChannel, futures::Vector{Future};
                         secs::Float64=0.1)

    n_alive = length(futures)

    @debug("Checking worker status.")

    for (i, f) in enumerate(futures)

        if timedwait(()->isready(f), secs) == :ok

            @warn("Worker $(i) seems to have finished prematurely.",
                  fetch(f))

            n_alive -= 1

        end

    end

    @debug("Workers online: $(n_alive)")

    # Only closes all queues if there are tasks to be completed
    if n_alive == 0 && isopen(tqueue)

        @warn("No live worker found. Will close queues and finish.")

        close(bqueue)

        close(tqueue)

        close(squeue)

    end

end
