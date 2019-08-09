using Cxx
using Libdl

push!(DL_LOAD_PATH, "/usr/local/lib")
push!(DL_LOAD_PATH, "./")

function ransac_example()

    for header in ["/usr/include/eigen3/",
                   "/usr/local/include/theia/libraries/vlfeat",
                   "/usr/include/suitesparse",
                   "/usr/local/include/theia/libraries/",
                   "/usr/local/include/theia/libraries/statx",
                   "/usr/local/include/theia/libraries/optimo",
                   "/usr/local/include",
                   "/usr/local/include/theia",
                   "/usr/include", "./"]
    
        addHeaderDir(header, kind=C_User)

    end

    Libdl.dlopen("libvlfeat.so", Libdl.RTLD_GLOBAL)
    Libdl.dlopen("libflann_cpp.so", Libdl.RTLD_GLOBAL)
    Libdl.dlopen("libakaze.so", Libdl.RTLD_GLOBAL)
    Libdl.dlopen("libceres.so", Libdl.RTLD_GLOBAL)
    Libdl.dlopen("libtheia.so", Libdl.RTLD_GLOBAL)
    Libdl.dlopen("libransac.so", Libdl.RTLD_GLOBAL)    

    cxx"""
run_ransac();

    """
    
end
