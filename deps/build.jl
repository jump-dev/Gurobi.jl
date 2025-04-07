# Copyright (c) 2015 Dahua Lin, Miles Lubin, Joey Huchette, Iain Dunning, and
#   contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

using Libdl

const DEPS_FILE = joinpath(@__DIR__, "deps.jl")

if isfile(DEPS_FILE)
    rm(DEPS_FILE)
end

if Int === Int32
    error(
        "Gurobi.jl does not support 32-bit Julia. Please install a 64-bit Julia.",
    )
end

function write_depsfile(path)
    open(DEPS_FILE, "w") do io
        println(io, "const libgurobi = \"$(escape_string(path))\"")
        return
    end
    return
end

const ALIASES =
    ["gurobi120", "gurobi110", "gurobi100", "gurobi95", "gurobi91", "gurobi90"]

function _try_local_install()
    paths_to_try = copy(ALIASES)
    for a in ALIASES
        root = get(ENV, "GUROBI_HOME", nothing)
        if root !== nothing
            if Sys.isunix()
                push!(paths_to_try, joinpath(root, "lib", "lib$a.so"),)
            end
            if Sys.iswindows()
                push!(paths_to_try, joinpath(root, "bin", "$a.$(Libdl.dlext)"))
            end
            if Sys.isapple()
                push!(paths_to_try, joinpath(root, "lib", "lib$a.dylib"))
            end
        end
        if Sys.isapple()  # gurobi uses .so on OS X for some reason
            push!(paths_to_try, string("lib$a.so"))
            push!(paths_to_try, string("lib$a.dylib"))
        end
    end
    for l in paths_to_try
        if Libdl.dlopen_e(l) != C_NULL
            write_depsfile(l)
            return true
        end
    end
    return false
end

function _print_GUROBI_HOME_help()
    version = "1000"
    println("""
    You should set the `GUROBI_HOME` environment variable to point to the
    install location then try again. For example (updating the path to the
    correct location if needed):
    ```
    # On Windows, this might be
    ENV["GUROBI_HOME"] = "C:\\\\Program Files\\\\gurobi$(version)\\\\win64\\\\"
    import Pkg
    Pkg.add("Gurobi")
    Pkg.build("Gurobi")

    # On OSX, this might be
    ENV["GUROBI_HOME"] = "/Library/gurobi$(version)/mac64/"
    import Pkg
    Pkg.add("Gurobi")
    Pkg.build("Gurobi")

    # On Unix, this might be
    ENV["GUROBI_HOME"] = "/opt/gurobi$(version)/linux64/"
    import Pkg
    Pkg.add("Gurobi")
    Pkg.build("Gurobi")
    ```
    **Note: your path may differ. Check which folder you installed the Gurobi
    binary in, and update the path accordingly.**
    """)
end

function diagnose_gurobi_install()
    println("""

    **Unable to locate Gurobi installation. Running some common diagnostics.**

    Gurobi.jl only supports the following versions:
    """)
    println.(" - ", ALIASES)
    println("""

    Did you download and install one of these versions from gurobi.com?
    Installing Gurobi.jl via the Julia package manager is _not_ sufficient!
    """)
    if haskey(ENV, "GUROBI_HOME")
        dir = joinpath(ENV["GUROBI_HOME"], Sys.isunix() ? "lib" : "bin")
        println("""
        Found GUROBI_HOME =  $(ENV["GUROBI_HOME"])

        Does this point to the correct install location?

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
            `libgurobiXXX.so`, `libgurobiXXX.dylib`, or `gurobiXXX.dll`.\n\n""")

            _print_GUROBI_HOME_help()
        catch ex
            if typeof(ex) <: SystemError
                println("""
                Aha! We tried looking in `$(dir)`, but something went wrong. Are
                you sure that your GUROBI_HOME environment variable is correct?
                When combined with the appropriate suffix (e.g., `lib` or
                `bin`, it needs to point to a valid directory.\n\n""")
                _print_GUROBI_HOME_help()
            else
                rethrow(ex)
            end
        end
    else
        try
            # Try to call `gurobi_cl`. This should work if Gurobi is on the
            # system path. If it succeeds, it will print out the version.
            io = IOBuffer()
            run(pipeline(`gurobi_cl --version`; stdout = io))
            seekstart(io)
            println("""
            We couldn't find the `GUROBI_HOME` environment variable, but we
            found this version of Gurobi on your `PATH`.

            $(read(io, String))
            Is this version one of the supported versions listed above? If so,
            we found the executable, but not the libraries we need. Follow the
            advice below to set the `GUROBI_HOME` environment variable. If not,
            you should edit your `PATH` to point to the correct version, or set
            the `GUROBI_HOME` environment variable.\n""")
            _print_GUROBI_HOME_help()
        catch
            println("""

            We could not find a version of Gurobi in your `PATH`, and we could
            not find the environment variable `GUROBI_HOME`.\n\n""")
            _print_GUROBI_HOME_help()
        end
    end
    return error(
        """
        Unable to locate Gurobi installation. If the advice above did not help,
        open an issue at https://github.com/jump-dev/Gurobi.jl and post the full
        print-out of this diagnostic attempt.
        """,
    )
end

if haskey(ENV, "GUROBI_JL_SKIP_LIB_CHECK")
    # We write a fake depsfile so Gurobi.jl is loadable but not usable.
    write_depsfile("__skipped_installation__")
elseif get(ENV, "JULIA_REGISTRYCI_AUTOMERGE", "false") == "true"
    # We write a fake depsfile so Gurobi.jl is loadable but not usable.
    write_depsfile("__skipped_installation__")
elseif get(ENV, "GUROBI_JL_USE_GUROBI_JLL", "true") == "false"
    # The user has asked to avoid Gurobi_jll
    found = _try_local_install()
    if !found
        diagnose_gurobi_install()
    end
else
    # We're using the artifact
    open(DEPS_FILE, "w") do io
        println(io, "# No libgurobi constant; we're using the Artifact.")
    end
end
