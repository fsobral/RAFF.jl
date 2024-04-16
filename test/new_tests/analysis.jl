using Plots, DelimitedFiles

function visualize!(file1::String,file2::String;legendas = ["metodo 1", "método 2"],param = "median_PT",saveas="new.pdf")
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
                push!(a,A[k,param_ind]) 
                push!(b,B[k,param_ind])
        end       
    end
    plot()
    if param=="niter"
        xlabel!("Problems")
        ylabel!("Number of iterations")
        bar!(a, label = "$(legendas[1])",alpha=0.5,fillcolor="purple")
        bar!(b, label = "$(legendas[2])",alpha=1.0,fillcolor="purple")
        savefig("$(saveas).pdf")
    else

        plot([1:length(a)],[a,b], label=["$(legendas[1])" "$(legendas[2])"], linewidth=2)
    end 
end

function analyse!(file1::String,file2::String;params = ["median_PT","memory"])
    A = readdlm(file1,',')
    B = readdlm(file2,',')

    conv_ind = findall(x->x=="convergence",A[1,:])[1]
    println(" 😏 Analise: ")
    println(" Os resultados de $(file1) indicam que o método convergiu em $(length(findall(A[2:end,conv_ind]))) problemas")
    println(" Os resultados de $(file2) indicam que o método convergiu em $(length(findall(B[2:end,conv_ind]))) problemas")

    error_ind = findall(x->x=="error",A[1,:])[1]
    ia = (findall(x->abs(x)<1.0e-3,A[2:end,error_ind])) .+ 1    
    ib = (findall(x->abs(x)<1.0e-3,B[2:end,error_ind])) .+ 1
    println(" Os resultados de $(file1) indicam que o método encontrou resposta com erro < 1.0e-3 em $(length(ia)) problemas.") 
    println(" Destes, o pior valor de erro encontrado foi $(maximum(A[ia,error_ind])) problemas, a média foi $(sum(A[ia,error_ind])/length(ia)) e o melhor foi $(minimum(A[ia,error_ind])).") 

    println(" Os resultados de $(file2) indicam que o método encontrou resposta com erro < 1.0e-3 em $(length(ib)) problemas.")
    println(" Destes, o pior valor de erro encontrado foi $(maximum(B[ib,error_ind])) problemas, a média foi $(sum(B[ib,error_ind])/length(ib)) e o melhor foi $(minimum(B[ib,error_ind])).") 


    for prs ∈ params 
        param_ind = findall(x->x== prs,A[1,:])[1]
        #println(param_ind)
        n = length(A[:,1])
        a = 0 
        b = 0 
        mp1 = 0.0
        mp2 = 0.0 
        nconv = 0
        for k=2:n
            if A[k,conv_ind]==true && B[k,conv_ind]==true
                nconv += 1
                mp1 = mp1+A[k,param_ind]
                mp2 = mp2+B[k,param_ind]
                if A[k,param_ind]<B[k,param_ind]
                    a += 1 
                else
                    b += 1 
                end
            end       
        end
        mp1 = mp1/nconv
        mp2 = mp2/nconv 
        println(" Considerando problemas em que os dois métodos convergeriam, o método dado em $(file1) obteve menor $(prs) em $(a) problemas.")
        println(" Considerando problemas em que os dois métodos convergeriam, o método dado em $(file2) obteve menor $(prs) em $(b) problemas.")
        println(" A média do indicador $(prs) no método dado em $(file1) foi de $(mp1).")
        println(" A média do indicador $(prs) no método dado em $(file2) foi de $(mp2).")
    end
    println(" 🏁 " )
e
