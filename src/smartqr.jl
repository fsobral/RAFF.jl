export cos_sin, create_vec, vector!, givens!, update!, solve_lower_vec!, solve_upper_vec!, solve_update_vec!

cos_sin(a::Float64, b::Float64)=
begin
#     if b == 0
#         c = 1
#         s = 0
#     elseif abs(b) > abs(a)
#         t = -(a / b)
#         s = 1 / sqrt(1 + (t ^ 2))
#         c = s * t
#     else
#         t = -(b / a)
#         c = 1 / sqrt(1 + (t ^ 2))
#         s = c * t
#     end
    c = a / (sqrt(a ^2 + b^2))
    s = - b / (sqrt(a ^2 + b^2))
    return c, s
end

create_vec(R)=
begin
    m, n = size(R)
    if n > m
        return []
    end
    z = Int64(n * (n + 1) / 2)
    vetor_R = zeros(z)
    return vetor_R
end

"""
    vector!(R, vetor_R) 

Overwrite the vector `vector_R` with the upper triangular matrix `R`. 

"""
vector!(R, vetor_R)=
begin
    m, n = size(R)
    if n > m
        return []
    end
    h = 1
    for i = 1 : n
        s = h + n - i
        vetor_R[h : s] = R[i, i : n]
        h = s + 1
    end
    return vetor_R
end

givens!(A::Array{Float64,1}, aux, aux2)=
begin    
    
    c, s = cos_sin(A[1], aux[1])
    aux2 .= A
    A .= A .* c .- s .* aux
    aux .= aux2 .* s .+ c .* aux
  
    return A, aux
end

"""
    update!(R::Array{Float64,1}, λ)

Receives the vector `R` that contains the upper triangular matrix of the qr factorization and a number `λ`. It computes a
new qr factorization of the upper triangular matrix with the square root of λ added in the main diagonal. It uses Givens rotations.

"""
update!(R::Array{Float64,1}, λ)=
begin
    z = length(R)   
    n = Int64((sqrt(1 + 8 * z) - 1) / 2)
    s_1 = 0
    aux = zeros(n)
    aux2 = zeros(n)
    
    for k = 1 : n
        aux[k] = sqrt(λ)
        s_1 = s_1 + n - k + 1 
        h_1 = s_1 - n + k 
        s = s_1
        h = h_1
        for t = k : n
            R[h : s], aux[t : n] = givens!(R[h : s], aux[t : n], aux2[t : n]) 
            h = s + 1 
            s = s + n - t 
        end
       
    end
    return R # = R_λ 
end

solve_lower_vec!(R_λ::Array{Float64,1}, b::Array{Float64,1}, z::Int64, n::Int64)=
begin
    s = 1
    for j = 1 : n 
        b[j] = b[j] / R_λ[s]
        s = s + 1
        for i = j + 1 : n
            b[i] = b[i] - R_λ[s] * b[j]
            s = s + 1
        end
    end
    return b
end

solve_upper_vec!(R_λ::Array{Float64,1}, b::Array{Float64,1}, z::Int64, n::Int64)=
begin
    s = z
    for k = n : -1 : 1
        sum = 0
        for i = n : -1 : k + 1
            sum = sum + b[i] * R_λ[s]
            s = s - 1 
        end
        b[k] = (b[k] - sum) / R_λ[s]
        s = s - 1
    end
    return b
end

"""
    solve_update_vec!(R_λ::Array{Float64,1}, b::Array{Float64,1})

Solves the system (R^T * R + λI)x = b, where R is given by the vector `R_λ` that contains the updated upper triangular
matrix R given by the qr factorization, and `b` that contains the right side of the system. `b` will be overwritten by the 
solution of the system.

"""
solve_update_vec!(R_λ::Array{Float64,1}, b::Array{Float64,1})= 
begin
    z = length(R_λ)
    n = length(b)
    
    b = solve_lower_vec!(R_λ, b, z, n)
    
    b = solve_upper_vec!(R_λ, b, z, n)
    
    return b
end
