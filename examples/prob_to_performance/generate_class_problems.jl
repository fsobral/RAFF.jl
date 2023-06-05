function generate!(probtype,nprob)
    if probtype == "parabola"
        for i=1:nprob
            l = rand([1.0:4.0;])
            a = rand([-10.0:0.5:10.0;])
            b = rand([-10.0:0.5:10.0;])
            c = rand([-10.0:0.5:10.0;])
            npts = rand([100:50:1000;])
            nout = round(Int,((rand([5:1:25;]))/100.0)*npts)

            build_problem("parabola",[-l,l],[a,b,c,npts,nout])
        end
    end
    if probtype == "cubic"
        for i=1:nprob
            l = rand([1.0:4.0;])
            a = rand([-10.0:0.5:10.0;])
            b = rand([-10.0:0.5:10.0;])
            c = rand([-10.0:0.5:10.0;])
            d = rand([-10.0:0.5:10.0;])
            npts = rand([100:50:1000;])
            nout = round(Int,((rand([5:1:25;]))/100.0)*npts)
            build_problem("cubic",[-l,l],[a,b,c,d,npts,nout])
        end
    end
    if probtype == "gaussian"
        for i=1:nprob
            l = rand([1.0:4.0;])
            a = rand([-10.0:0.5:10.0;])
            b = rand([-10.0:0.5:10.0;])
            c = rand([-10.0:0.5:10.0;])
            npts = rand([100:50:1000;])
            nout = round(Int,((rand([5:1:25;]))/100.0)*npts)
            build_problem("gaussian",[-l,l],[a,b,c,npts,nout])
        end
    end
    if probtype == "log"
        for i=1:nprob
            l = rand([2.0:5.0;])
            a = rand([1.0:0.5:10.0;])
            b = rand([1.0:0.5:10.0;])
            c = rand([1.0:0.5:10.0;])
            npts = rand([100:50:1000;])
            nout = round(Int,((rand([5:1:25;]))/100.0)*npts)
            build_problem("log",[1,l],[a,b,c,npts,nout])
        end
    end
     if probtype == "trig"
        for i=1:nprob
            l = rand([2.0:5.0;])
            a = rand([1.0:0.5:10.0;])
            b = rand([1.0:0.5:10.0;])
            c = rand([1.0:0.5:10.0;])
            npts = rand([100:50:1000;])
            nout = round(Int,((rand([5:1:25;]))/100.0)*npts)
            build_problem("trig",[1,l],[a,b,c,npts,nout])
        end
    end
    



end

generate_all() = begin
    #generate!("parabola",250)
    #generate!("gaussian",250)
    #generate!("cubic",250)
    generate!("trig",250)
end
