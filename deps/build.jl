depsfile = joinpath(dirname(@__FILE__),"deps.jl")
if isfile(depsfile)
    rm(depsfile)
end

function write_depsfile(path)
    f = open(depsfile,"w")
    println(f,"const libgurobi = \"$path\"")
    close(f)
end

aliases = ["gurobi60","gurobi56","gurobi55","gurobi51"]

paths_to_try = [aliases]

for a in aliases
    if haskey(ENV, "GUROBI_HOME")
        @unix_only push!(paths_to_try, joinpath(ENV["GUROBI_HOME"],"lib",string("lib",a,".so")))
        @windows_only push!(paths_to_try, joinpath(ENV["GUROBI_HOME"],"bin",string(a,".",Sys.dlext)))
    end
    # gurobi uses .so on OS X for some reason
    @osx_only push!(paths_to_try, string("lib$a.so"))
end

for l in paths_to_try
    d = dlopen_e(l)
    if d != C_NULL
        write_depsfile(l)
        break
    end
end
