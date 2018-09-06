#robust algebric fitting function (RAFF)

#dependencies
using ForwardDiff

function LMlovo(model::Function,data::Array{Float64,2},n::Int,p::Int)
    # model = the input model 
    # data =  the dataset 
    # n = number of variables in the model
    # p = trusted points

    #defining closures for derivative and initializations
    model_cl(x)=model(x,t)
    grad_model_cl(x)=ForwardDiff.gradient(model_cl, x)
    t=0.0
    npun=length(data[:,1])
    
    #All evaluated instance depends of a sort function 
    SortFun(V::Vector)=begin
        aux=0
        vaux=0.0
        ind=[1:1:npun;]
        for i=1:p
            for j=i+1:npun
                if (V[i]>V[j])
                    aux=ind[j]
                    ind[j]=ind[i]
                    ind[i]=aux
                    vaux=V[j]
                    V[j]=V[i]
                    V[i]=vaux
                end
            end
        end
        return ind[1:1:p]
    end
    
    #main functions
    LovoFun(x)=begin #return a ordered set index and lovo value 
        F=zeros(npun)
        for i=1:npun 
            F[i]=(model(x,data[i,1])-data[i,2])^2
        end
        return SortFun(F),sum(F[SortFun(F)])
    end
    ResFun(x,ind)=begin #return the residue and Jacobian of residue 
        r=zeros(p)
        rJ=zeros(p,n)
        k=1
        for i in ind
            r[k]=model(x,data[i,1])-data[i,2]
            t=data[i,1]
            rJ[k,:]=grad_model_cl(x)
            k=k+1
        end
        return r,rJ 
    end
    
    #Levenberg-Maquardt algorithm
    I=eye(n)
    ε=10.0^(-4)
    λ_up=2.0
    λ_down=2.0
    λ=1.0
    x=zeros(n) #initial point
    d=zeros(n)
    y=zeros(n)
    (ind_lovo,val_lovo)=LovoFun(x)
    (val_res,jac_res)=ResFun(x,ind_lovo)
    grad_lovo=jac_res'*val_res
    safecount=1
    while norm(grad_lovo,2)>=ε && safecount<200
        G=jac_res'*jac_res+λ*I 
        F=qrfact(G)     
        ad=try
            y=F[:Q]\(-grad_lovo)
            d=F[:R]\y
        catch
            "error"
        end
        if ad=="error" #restarting if lapack fails
            d=-grad_lovo 
            x=rand(n)
        else 
            d=ad
        end
        #d=G\(-grad_lovo)
        xnew=x+d
        (ind_lovo_new,val_lovo_new)=LovoFun(xnew)
        if  val_lovo_new<=val_lovo 
            x=copy(xnew)
            val_lovo=copy(val_lovo_new)
            ind_lovo=copy(ind_lovo_new)
            λ=λ/(λ_down)
            (val_res,jac_res)=ResFun(x,ind_lovo)
            grad_lovo=jac_res'*val_res
        else
            λ=λ*(λ_up)
        end
        safecount+=1
    end
    if safecount==200
    #    println("no solution was founded in $safecount iterations")
        return 0,x,safecount,p
    else
    #    println("solution founded::   $x " )
    #    println("number of iterations:: $(safecount)")
        return 1,x,safecount,p
    end
end


function raff(model,data,n)
    pliminf=Int(round(length(data[:,1])/2.0))
    plimsup=length(data[:,1])
    v=Array{Any,1}(plimsup-pliminf+1)
    k=1
    for i=pliminf:plimsup
        v[k]=LMlovo(model,data,n,i)
        k+=1
    end
    lv=length(v)
    votsis=zeros(lv)
    for i=1:lv
        for j=1:lv
            if norm(v[i][2]-v[j][2])<10.0^(-3)
                votsis[i]+=1
            end
        end
    end
    mainind=findlast(x->x==maximum(votsis),votsis)
    display(v[mainind])
end


##############################################
###########testing############################
##############################################

#adjustment model
#function model(x,t)
#    return x[1]*exp(t*x[2])
#end

#exact answer 2exp(-0.5t)

#data=[-1.0 3.2974425414002564;
#     -0.75 2.9099828292364025;
#     -0.5 2.568050833375483;
#     -0.25 2.2662969061336526;
#     0.0 2.0;
#     0.25 1.764993805169191;
#     0.5 1.5576015661428098;
#     0.75 1.5745785575819442; #noise
#     1.0 1.2130613194252668;
#     1.25 1.0705228570379806;
#     1.5 0.9447331054820294;
#     1.75 0.8337240393570168;
#     2.0 0.7357588823428847;
#     2.25 0.6493049347166995;
#     2.5 0.5730095937203802;
#     2.75 0.5056791916094929;
#     3.0 0.44626032029685964;
#     3.25 0.5938233504083881; #noise 
#     3.5 0.3475478869008902;
#     3.75 0.30670993368985694;
#     4.0 0.5706705664732254;] #noise
#
##LMlovo(model,data,2,18)      

#raff(model,data,2)