@testset "Smart qr" begin

@testset "Test vector" begin
   
    M = [1 2; 4 5]
    v = create_vec(M)
    v = vector!(M, v)
    
    @test length(v) == 3
    @test v[1] == M[1,1]
    @test v[2] == M[1,2]
    @test v[3] == M[2,2]

    # Test a matriz which m > n
    
    M = [1 2; 4 5; 6 7]
    v = create_vec(M)
    v = vector!(M, v)
    
    @test length(v) == 3
    @test v[1] == M[1,1]
    @test v[2] == M[1,2]
    @test v[3] == M[2,2]

    # Test an invalid case
    
    M = [1 2 3; 4 5 6]
    
    @test length(create_vec(M)) == 0
        
end

@testset "Test givens" begin
   
    a = [1.0, 0, 0, 0]
    aorig = copy(a)
    b = [1.0, 0, 0, 0]
    aux2 = zeros(length(a))
    
    givens!(a, b, aux2)
    
    @test( b == zeros(length(b)))
    @test( aux2 == aorig)

    a = [1.0, 0, 0, 0]
    b = [-1.0, 0, 0, 0]
    aux2 = zeros(length(a))
    
    givens!(a, b, aux2)
    
    @test( b == zeros(length(b)))
    
    a = [1.0, 1.0, 1.0, 1.0]
    b = [-1.0, 0, 0, 0]
    aux2 = zeros(length(a))
    
    givens!(a, b, aux2)
    
    @test(b[2] == b[3] == b[4])
    # Fazer as contas
    

end

@testset "Test update" begin
   
    R = [1.0 2 3 4 5; 0 6 7 8 9; 0 0 10 11 12; 0 0 0 13 14; 0 0 0 0 15]
    
    λ = 4.0
    
     jQ, jR = qr([R; sqrt(λ) * Matrix(I(5))])
    #jR, = QR_fac([R; sqrt(λ) * Matrix(I(5))])
    w = create_vec(jR)
    w = vector!(jR, w)
    
    v = create_vec(R)
    v = vector!(R, v)
    
    update!(v, λ)
    
    @test(norm(v) ≈ norm(w))
    
    for i = 1 : length(v)
        @test(abs(v[i]) ≈ abs(w[i]))
    end
    
end

@testset "Test update" begin
   
    R = [sqrt(3) -17.678 0 485 5.2 14; 0 694 -10 π 197 94; 0 0 -5 11 -24 1000; 0 0 0 19 2 87; 0 0 0 0 5 0; 0 0 0 0 0 -344]
    
    λ = 4.0

    #função de julia
    jQ, jR = qr([R; sqrt(λ) * Matrix(I(6))])
    w = create_vec(jR)
    w = vector!(jR, w)
    
    #minha função
    v = create_vec(R)
    v = vector!(R, v)
    
    update!(v, λ)
    
    @test(norm(v) ≈ norm(w))
    
    for i = 1 : length(v)
        @test(abs(v[i]) ≈ abs(w[i]))
    end
    
end

@testset "Test solution" begin
   
    R = [1.0 2 3 4 5; 0 6 7 8 9; 0 0 10 11 12; 0 0 0 13 14; 0 0 0 0 15]
    
    λ = 4.0
    
    b = [0, 5.0, -7, 2, 3]
    
    #função de julia
    jQ, jR = qr([R; sqrt(λ) * Matrix(I(5))])
    
    c = (transpose(jR) * jR) \ b
    
    
    #minha função
    v = create_vec(R)
    v = vector!(R, v)
    
    update!(v, λ)
    
    solve_update_vec!(v, b)
    
    for i = 1 : length(b)
        @test(b[i] ≈ c[i])
    end
    
    
end

@testset "Test solution" begin
   
    R = [sqrt(3) -17.678 0 485 5.2 14; 0 694 -10 π 197 94; 0 0 -5 11 -24 1000; 0 0 0 19 2 87; 0 0 0 0 5 0; 0 0 0 0 0 -344]
    
    λ = 4.0
    
    b = [123, 9.0, -0.15, 49, 3, sqrt(14)]
    
    #função de julia
    jQ, jR = qr([R; sqrt(λ) * Matrix(I(6))])
    
    c = (transpose(jR) * jR) \ b
    
    #minha função
    v = create_vec(R)
    v = vector!(R, v)
    
    update!(v, λ)
    
    solve_update_vec!(v, b)
    
    for i = 1 : length(b)
        @test(b[i] ≈ c[i])
    end
    
    
    
end
end
