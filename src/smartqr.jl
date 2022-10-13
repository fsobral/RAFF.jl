export cs, create_vec, vector!, givens3, update3, solve_lower_vec, solve_upper_vec, solve_update_vec

cs(a, b)=
begin
#     if b ==0
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

#transforma uma matriz triangular superior em um vetor
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

givens3(A, aux, aux2)=
begin    
    
    #zeramos o elemento 1 da linha recebida
    c, s = cs(A[1], aux[1])
    aux2 .= A
    A .= A .* c .- s .* aux
    aux .= aux2 .* s .+ c .* aux
  
    return A, aux
end

#sobreescrevemos a fatoração QR recebida
update3(R, λ)=
begin
    z = length(R)   
    n = Int64((sqrt(1 + 8 * z) - 1) / 2)
    s_1 = 0
    aux = zeros(n)
    aux2 = zeros(n)
    
    for k = 1 : n
        aux[k] = sqrt(λ)
        s_1 = s_1 + n - k + 1 #ultimo indice da linha k 
        h_1 = s_1 - n + k #primeiro indice da linha k
        s = s_1
        h = h_1
        for t = k : n
            R[h : s], aux[t : n] = givens3(R[h : s], aux[t : n], aux2[t : n]) #vai zerar o primiero elemento de aux
            h = s + 1 #vai para a proxima linha
            s = s + n - t #vai para o ultimo elemento da proxima linha
        end
       
    end
    return R # = R_λ 
end

solve_lower_vec(R_λ, b, z, n)=
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

solve_upper_vec(R_λ, b, z, n)=
begin
    s = z
    for k = n : -1 : 1
        soma = 0
        for i = n : -1 : k + 1
            soma = soma + b[i] * R_λ[s]
            s = s - 1 
        end
        b[k] = (b[k] - soma) / R_λ[s]
        s = s - 1
    end
    return b
end

solve_update_vec(R_λ, b)= #R ja alterada
begin
    z = length(R_λ)
    n = length(b)
    
    b = solve_lower_vec(R_λ, b, z, n)
    
    b = solve_upper_vec(R_λ, b, z, n)
    
    return b
end
