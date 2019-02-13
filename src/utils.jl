"""

    eliminate_local_min!(sols::Vector{RAFFOutput})

Check if the function value of the solution found by smaller values of
`p` is not greater when compared with larger ones. This certainly
indicates that a local minimizer was found by the smaller `p`.

"""
function eliminate_local_min!(model::Function, data::Array{Float64, 2},
                              sols::Vector{RAFFOutput})

    lv = length(sols)

    sec_sol_ind = -1

    # Start from the largest p
    for i = lv:-1:1

        (sols[i].status != 1) && continue
        
        for j = i - 1:-1:1

            if (sols[j].status == 1) && (sols[j].f > sols[i].f)

                @debug("  Possible local minimizer for p = $(sols[j].p) " *
                       "with f = $(sols[j].f). Removing it.")
                
                sols[j] = RAFFOutput(0)

            end

        end

    end

    # Test the maximum 'p'

    nump, = size(data)
    
    maxp = sols[lv]
        
    if (lv > 1) && (maxp.status == 1)

        i = lv - 1

        # while (i > 0) && (sols[i].status != 1)

        #     i -= 1

        # end

        bestf = maxp.f
        j = 0
        
        while (i > 0)

            if  (sols[i].status == 1) && (sols[i].f < bestf)

                bestf = sols[i].f

                j = i

            end

            i -= 1

        end

        if j == 0
            
            @debug("  Nobody to reject p = $(maxp.p).")

        else

            sec_sol = sols[j]
            
            @debug("  p = $(sec_sol.p) will try to reject p = $(maxp.p).")

            nmin = 0

            for i = 1:nump

                y  = data[i, end]
                t  = @view data[i, 1:(end - 1)]
                y1 = model(maxp.solution, t)
                y2 = model(sec_sol.solution, t)
                
                (abs(y - y1) < abs(y - y2)) && (nmin += 1)

            end

            if nmin < lv / 2

                @debug("  nmin = $(nmin). Rejecting p = $(maxp.p).")

                sols[lv] = RAFFOutput(0)

            end

        end
        
    end
    
end


"""

This function is an auxiliary function. It finds the `p` smallest
values of vector `V` and brings them to the first `p` positions. The
indexes associated with the `p` smallest values are stored in `ind`.

"""
function SortFun!(V::Vector{Float64}, ind::Vector{Int}, p::Int)

    # If p is invalid, returns an empty view
    (p <= 0) && (return @view(ind[1:p]), @view(V[1:p]))
    
    npun = length(ind)
    
    for i = 1:npun
        ind[i] = i
    end
    
    for i = 1:p
        
        for j = i + 1:npun
            
            if V[i] > V[j]
                
                aux = ind[j]
                ind[j] = ind[i]
                ind[i] = aux
                
                vaux = V[j]
                V[j] = V[i]
                V[i] = vaux
            end
        end
    end

    return @view(ind[1:p]), @view(V[1:p])
    
end

function setRaffOutputLevel(level::LogLevel)

    global raff_logger = ConsoleLogger(stdout, level)

end

function setLMOutputLevel(level::LogLevel)

    global lm_logger = ConsoleLogger(stdout, level)

end
