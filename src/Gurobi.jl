module Gurobi

import Pkg

const _DEPS_FILE = joinpath(dirname(@__FILE__), "..", "deps", "deps.jl")

if isfile(_DEPS_FILE)
    include(_DEPS_FILE)
elseif Sys.islinux()
    # If there is no _DEPS_FILE and we're on linux, use the Artifact
    # installation.
    const libgurobi = joinpath(
        Pkg.Artifacts.artifact"gurobilinux64",
        "gurobi950/linux64/lib/libgurobi95.so",
    )
else
    error("""
        Gurobi not properly installed. Please run Pkg.build(\"Gurobi\"). For
        more information go to https://github.com/jump-dev/Gurobi.jl
    """)
end

using CEnum

const _GUROBI_VERSION = if libgurobi == "julia_registryci_automerge"
    VersionNumber(9, 5, 0)
else
    let
        majorP, minorP, technicalP = Ref{Cint}(), Ref{Cint}(), Ref{Cint}()
        ccall(
            (:GRBversion, libgurobi),
            Cvoid,
            (Ptr{Cint}, Ptr{Cint}, Ptr{Cint}),
            majorP, minorP, technicalP,
        )
        VersionNumber(majorP[], minorP[], technicalP[])
    end
end

function _is_patch(x::VersionNumber, reference::VersionNumber)
    return x.major == reference.major && x.minor == reference.minor
end

if _is_patch(_GUROBI_VERSION, v"9.0")
    include("gen90/ctypes.jl")
    include("gen90/libgrb_common.jl")
    include("gen90/libgrb_api.jl")
elseif _is_patch(_GUROBI_VERSION, v"9.1")
    include("gen91/ctypes.jl")
    include("gen91/libgrb_common.jl")
    include("gen91/libgrb_api.jl")
elseif _is_patch(_GUROBI_VERSION, v"9.5")
    include("gen95/ctypes.jl")
    include("gen95/libgrb_common.jl")
    include("gen95/libgrb_api.jl")	
else
    error("""
    You have installed version $_GUROBI_VERSION of Gurobi, which is not
    supported by Gurobi.jl. We require Gurobi version 9.0 or 9.1 or 9.5. 

    After installing a supported version of Gurobi, run:

        import Pkg
        Pkg.rm("Gurobi")
        Pkg.add("Gurobi")

    Make sure you set the environment variable `GUROBI_HOME` following
    the instructions in the Gurobi.jl README, which is available at
    https://github.com/jump-dev/Gurobi.jl.

    If you have a newer version of Gurobi installed, changes may need to be made
    to the Julia code. Please open an issue at
    https://github.com/jump-dev/Gurobi.jl.
    """)
end

include("MOI_wrapper/MOI_wrapper.jl")
include("MOI_wrapper/MOI_callbacks.jl")
include("MOI_wrapper/MOI_multi_objective.jl")
include("MOI_wrapper/MOI_indicator_constraint.jl")

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
