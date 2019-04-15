# This script generates random test problems

A = ["(x,t)->x[1]*exp(t*x[2])"     2;
     "(x,t)->x[1]*t+x[2]"          2;
     "(x,t)->x[1]*t^2+x[2]*t+x[3]" 3]

number_of_examples = 2

number_of_functions = length(A[:,1])

tmin = - 10.0

tmax = 10.0

v_points = [10:10:1000;]

list = open("list.dat", "w")

for i = 1:number_of_functions
    
    n = A[i,2]
    
    for j = 1:number_of_examples
        
        np = rand(v_points)
        
        p = rand([Int(round(np / 2)):1:np;])
        
        modelStr = A[i, 1]
        
        model = eval(parse(modelStr))

        solFilename = "sol_example_$(i)_$(j)_$(np)_$(p).dat"

        datFilename = "example_$(i)_$(j)_$(np)_$(p).dat"
                
        println(list, datFilename, "  ", solFilename)

        generate_test_problems(datFilename, solFilename, model,
                             modelStr, n, np, p, tmin = tmin,
                             tmax = tmax)

    end
    
end

close(list)

