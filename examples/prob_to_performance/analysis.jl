using Plots, DelimitedFiles

function visualize!(file1::String,file2::String;param = "median_PT")
    A = readdlm(file1,',')
    B = readdlm(file2,',')
    conv_ind = findall(x->x=="convergence",A[1,:])[1]
    println(conv_ind)
    error_ind = findall(x->x=="error",A[1,:])[1]
    println(error_ind)
    param_ind = findall(x->x== param,A[1,:])[1]
    println(param_ind)
    n = length(A[:,1])
    a = [] 
    b = [] 
    for k=2:n
        if A[k,conv_ind]==true && B[k,conv_ind]==true
            if max(A[k,error_ind],B[k,error_ind])<1.0e-3
                push!(a,A[k,param_ind]) 
                push!(b,B[k,param_ind])
            end
        end       
    end
    plot([1:length(a)],[a,b])
end

function analyse!(file1::String,file2::String;param = "median_PT")
    A = readdlm(file1,',')
    B = readdlm(file2,',')

    conv_ind = findall(x->x=="convergence",A[1,:])[1]
    println("↪ Os resultados de $(file1) indicam que o método convergiu em $(length(findall(A[2:end,conv_ind]))) problemas")
    println("↪ Os resultados de $(file2) indicam que o método convergiu em $(length(findall(B[2:end,conv_ind]))) problemas")

    error_ind = findall(x->x=="error",A[1,:])[1]
    println("↪ Os resultados de $(file1) indicam que o método encontrou resposta com erro < 1.0e-3 em $(length(findall(x->abs(x)<1.0e-3,A[2:end,error_ind]))) problemas") 
    println("↪ Os resultados de $(file2) indicam que o método encontrou resposta com erro < 1.0e-3 em $(length(findall(x->abs(x)<1.0e-3,B[2:end,error_ind]))) problemas")

    param_ind = findall(x->x== param,A[1,:])[1]
    println(param_ind)
    n = length(A[:,1])
    a = 0 
    b = 0 
    for k=2:n
        if A[k,conv_ind]==true && B[k,conv_ind]==true
            if A[k,error_ind]<1.0e-3
                a += 1 
            else
                b += 1 
            end
        end       
    end
    println("↪ Considerando problemas em que os dois métodos convergeriam, o método dado em $(file1) obteve menor $(param) em $(a) problemas")
    println("↪ Considerando problemas em que os dois métodos convergeriam, o método dado em $(file2) obteve menor $(param) em $(b) problemas")
end
