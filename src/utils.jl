"""

This function is an auxiliary function. It finds the `p` smallest
values of vector `V` and brings them to the first `p` positions. The
indexes associated with the `p` smallest values are stored in `ind`.

"""
function SortFun!(V::Vector{Float64}, ind::Vector{Int}, p::Int)
    
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
