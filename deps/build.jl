using Compat

@static if VERSION >= v"0.7"
    using Libdl
end

const DEPS_FILE = joinpath(@__DIR__, "deps.jl")

if isfile(DEPS_FILE)
    rm(DEPS_FILE)
end

function write_depsfile(path)
    if Compat.Sys.iswindows()
        # When `path` gets written out to a file, it will escape any
        # backslashes, so we need to doubly escape them. If your path uses
        # forward slashes, this operation won't do anything.
        path = replace(path, "\\" => "\\\\")
    end
    open(DEPS_FILE, "w") do io
        println(io, "const libgurobi = \"$(path)\"")
    end
end

const ALIASES = [
    "gurobi81", "gurobi80",
    "gurobi75", "gurobi70"
]

paths_to_try = copy(ALIASES)

for a in ALIASES
    if haskey(ENV, "GUROBI_HOME")
        if Compat.Sys.isunix()
            push!(paths_to_try, joinpath(ENV["GUROBI_HOME"], "lib", string("lib", a, ".so")))
        end
        if Compat.Sys.iswindows()
            push!(paths_to_try, joinpath(ENV["GUROBI_HOME"], "bin", string(a, ".", Libdl.dlext)))
        end
        if Compat.Sys.isapple()
            push!(paths_to_try, joinpath(ENV["GUROBI_HOME"], "lib", string("lib", a, ".dylib")))
        end
    end
    # gurobi uses .so on OS X for some reason
    if Compat.Sys.isapple()
        push!(paths_to_try, string("lib$a.so"))
        push!(paths_to_try, string("lib$a.dylib"))
    end
end

found = false
for l in paths_to_try
    d = Libdl.dlopen_e(l)
    if d != C_NULL
        global found = true
        write_depsfile(l)
        break
    end
end

if !found && !haskey(ENV, "GUROBI_JL_SKIP_LIB_CHECK")
    error("Unable to locate Gurobi installation. Note that this must be downloaded separately from gurobi.com")
end
