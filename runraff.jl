include("raff.jl")
N=readdlm("list.dat")
(m,n)=size(N)
for i=1:m 
    data=readdlm(N[i,1])[:,[1,2]]
    aux=readdlm(N[i,2])
    n=aux[1,1]
    tsol=""
    for k=1:n 
        tsol=tsol*(aux[2,k])#aux[2,] is substring
    end
    model=eval(parse(aux[3,1]))
    println(eval(parse(tsol)))
    raff(model,data,n)
end 
