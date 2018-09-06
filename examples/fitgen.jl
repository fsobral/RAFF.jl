A=["(x,t)->x[1]*exp(t*x[2])" 2;
"(x,t)->x[1]*t+x[2]" 2;
"(x,t)->x[1]*t^2+x[2]*t+x[3]" 3]

number_of_example=2
number_of_function=length(A[:,1])
tmin=-10.0
tmax=10.0
v_points=[10:10:1000;]
list=open("list.dat","w")
for i=1:number_of_function
    n=A[i,2]
    for j=1:number_of_example
        np=rand(v_points)
        p=rand([Int(round(np/2)):1:np;])
        #generate parameters x (solution)
        x=10.0*randn(n)
        sol=open("sol_example_$(i)_$(j)_$(np)_$(p).dat","w")
        println(sol,n) #number of variables
        println(sol,x) #parameters 
        println(sol,A[i,1]) #function expression
        close(sol)
        #
        #generate (ti,yi) where tmin<=t_i<=tmax (data)
        t=[tmin:(tmax-tmin)/(np-1):tmax;]
        y=zeros(length(t))
        data=open("example_$(i)_$(j)_$(np)_$(p).dat", "w")
        f=eval(parse(A[i,1]))
        v=rand([1:1:np;],np-p)
        for k=1:np
            y[k]=f(x,t[k])
            if k in v 
                noise=randn()*rand([1.0,2.0,3.0])
                println(data,t[k]," ",y[k]+noise,"     ",noise)
            else
                println(data,t[k]," ",y[k],"     ",0.0)
            end
        end
        close(data)
        
        println(list,"example_$(i)_$(j)_$(np)_$(p).dat", "  ","sol_example_$(i)_$(j)_$(np)_$(p).dat")
        #        
    end  
end
close(list)

