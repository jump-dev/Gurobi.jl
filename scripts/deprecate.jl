# This file was used to create the list of deprecations when moving from v0.8.1
# to v0.9.0.

using Gurobi

io = open("deprecated_functions.jl", "w")
print(io, """
const _DEPRECATED_ERROR_MESSAGE = \"\"\"
The C API of Gurobi.jl has been rewritten to expose the complete C API, and
all old functions have been removed.

For example:

    model = Gurobi.Optimizer()
    stat = Gurobi.get_status_code(model.inner)

is now:

    model = Gurobi.Optimizer()
    valueP = Ref{Cint}()
    ret = GRBgetintattr(model, "Status", valueP)
    if ret != 0
        # Do something because the call failed
    end
    stat = valueP[]

The new API is more verbose, but the names and function arguments are now
identical to the C API, documentation for which is available at:
https://www.gurobi.com/documentation/9.0/refman/c_api_details.html

To revert to the old API, use:

    import Pkg
    Pkg.add(Pkg.PackageSpec(name = \"Gurobi\", version = v\"0.8.1\"))

Then restart Julia for the change to take effect.
\"\"\"
""")

exported_names = Base.names(Gurobi; all = false)
for name in Base.names(Gurobi; all = true)
    foo = getfield(Gurobi, name)
    if !(foo isa Function)
        continue
    elseif any(startswith.(Ref(string(foo)), ["#", "@", "_"]))
        continue
    end
    println(io, "$(foo)(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)")
    if name in exported_names
        println(io, "export $(foo)")
    end
    println(io)
end
close(io)
