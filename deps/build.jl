using Libdl

const DEPS_FILE = joinpath(@__DIR__, "deps.jl")

if isfile(DEPS_FILE)
    rm(DEPS_FILE)
end

function write_depsfile(path)
    if Sys.iswindows()
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
        if Sys.isunix()
            push!(paths_to_try, joinpath(ENV["GUROBI_HOME"], "lib", string("lib", a, ".so")))
        end
        if Sys.iswindows()
            push!(paths_to_try, joinpath(ENV["GUROBI_HOME"], "bin", string(a, ".", Libdl.dlext)))
        end
        if Sys.isapple()
            push!(paths_to_try, joinpath(ENV["GUROBI_HOME"], "lib", string("lib", a, ".dylib")))
        end
    end
    # gurobi uses .so on OS X for some reason
    if Sys.isapple()
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

function diagnose_gurobi_install()
    println("""
    Unable to locate Gurobi installation. Running some common diagnostics.

    Gurobi.jl only supports the following versions:
    """)
    println.(" - ", ALIASES)
    println("""

    Did you download and install one of these versions from gurobi.com?

    """)
    if haskey(ENV, "GUROBI_HOME")
        dir = joinpath(ENV["GUROBI_HOME"], Sys.isunix() ? "lib" : "bin")
        println("""
        Found GUROBI_HOME =  $(ENV["GUROBI_HOME"])

        Does this point to the correct install location?
        - on Windows, this might be `C:\\Program Files\\gurobi810\\win64\\`
        - alternatively, on Windows, this might be `C:/Program Files/gurobi810/win64/`
        - on OSX, this might be `/Library/gurobi810/mac64/`
        - on Unix, this might be `/home/my_user/gurobi810/linux64/`

        Note: this has to be a full path, not a path relative to your current
        directory or your home directory.

        We're going to look for the Gurobi library in this directory:
            $(dir)

        That directory has the following files:
        """)
        try
            for file in readdir(dir)
                println(" - ", joinpath(dir, file))
            end
            println("""

            We were looking for (but could not find) a file named like
            `libgurobiXXX.so`, `libgurobiXXX.dylib`, or `gurobiXXX.dll`. You
            should update your GUROBI_HOME environment variable to point to the
            correct location.
            """)
        catch ex
            if typeof(ex) <: SystemError
                println("""
                Aha! We tried looking in `$(dir)`, but something went wrong. Are
                you sure that your GUROBI_HOME environment variable is correct?
                When combined with the appropriate suffix (e.g., `lib` or
                `bin`, it needs to point to a valid directory.
                """)
            else
                rethrow(ex)
            end
        end
    else
        try
            println("""
            Looking for a version of Gurobi in your path:
            """)
            # Try to call `gurobi_cl`. This should work if Gurobi is on the
            # system path. If it succeeds, it will print out the version.
            run(`gurobi_cl --version`)
            println("""

            We couldn't find the `GUROBI_HOME` environment variable, but we
            found this version of Gurobi on your path. Is this version one of
            the supported versions listed above? If not, you should edit your
            `PATH` to point to the correct version.
            """)
        catch
            println("""

            We could not find a version of Gurobi in your path, and we could
            not find the environment variable `GUROBI_HOME`. You should set 
            the `GUROBI_HOME` environment variable to point to the install
            location. For example:
            - on Windows, this might be `C:\\Program Files\\gurobi810\\win64\\`
            - on OSX, this might be `/Library/gurobi810/mac64/`
            - on Unix, this might be `/opt/gurobi810/linux64/`

            Alternatively, you can add the Gurobi install directory to your path.
            """)
        end
    end
end

if !found && !haskey(ENV, "GUROBI_JL_SKIP_LIB_CHECK")
    diagnose_gurobi_install()
    error("""
    Unable to locate Gurobi installation. If the advice above did not help,
    open an issue at https://github.com/JuliaOpt/Gurobi.jl and post the full
    print-out of this diagnostic attempt.
    """)
end
