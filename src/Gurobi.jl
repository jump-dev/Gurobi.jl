module Gurobi

const _DEPS_FILE = joinpath(dirname(@__FILE__), "..", "deps", "deps.jl")
if isfile(_DEPS_FILE)
    include(_DEPS_FILE)
else
    error("""
        Gurobi not properly installed. Please run Pkg.build(\"Gurobi\"). For
        more information go to https://github.com/jump-dev/Gurobi.jl
    """)
end

using CEnum

include("gen/ctypes.jl")
include("gen/libgrb_common.jl")
include("gen/libgrb_api.jl")

const _GUROBI_VERSION = VersionNumber(
    "$GRB_VERSION_MAJOR.$GRB_VERSION_MINOR.$(GRB_VERSION_TECHNICAL)"
)

if !(v"9.0.0" <= _GUROBI_VERSION < v"9.1")
    error(
        "You have installed version $_GUROBI_VERSION of Gurobi, which is not " *
        "supported by Gurobi.jl. If the version change was breaking, changes " *
        "may need to be made to the Julia code. Please open an issue at " *
        "https://github.com/jump-dev/Gurobi.jl."
    )
end

include("MOI_wrapper.jl")
include("MOI_callbacks.jl")
include("MOI_multi_objective.jl")

# Gurobi exports all `GRBXXX` symbols. If you don't want all of these symbols in
# your environment, then use `import Gurobi` instead of `using Gurobi`.

for sym in names(@__MODULE__, all=true)
    sym_string = string(sym)
    if startswith(sym_string, "GRB")
        @eval export $sym
    end
end

include("deprecated_functions.jl")

end
