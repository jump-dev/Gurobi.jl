# Copyright (c) 2015 Dahua Lin, Miles Lubin, Joey Huchette, Iain Dunning, and
#   contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

module Gurobi

import LazyArtifacts

# deps.jl file is always built via `Pkg.build`, even if we didn't find a local
# install and we want to use the artifact instead. This is so Gurobi.jl will be
# recompiled if we update the file. See issue #438 for more details.
include(joinpath(dirname(@__FILE__), "..", "deps", "deps.jl"))

if isdefined(@__MODULE__, :libgurobi)
    # deps.jl must define a local installation.
elseif Sys.islinux()
    # Let's use the artifact instead.
    const libgurobi = joinpath(
        LazyArtifacts.artifact"gurobilinux64",
        "gurobi1000/linux64/lib/libgurobi100.so",
    )
else
    error("""
        Gurobi not properly installed. Please run Pkg.build(\"Gurobi\"). For
        more information go to https://github.com/jump-dev/Gurobi.jl
    """)
end

const _GUROBI_VERSION = if libgurobi == "__skipped_installation__"
    # The deps file is fake, with the intention to make Gurobi.jl loadable but
    # not usable.
    VersionNumber(10, 0, 0)
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
    include("gen90/libgrb_common.jl")
    include("gen90/libgrb_api.jl")
elseif _is_patch(_GUROBI_VERSION, v"9.1")
    include("gen91/libgrb_common.jl")
    include("gen91/libgrb_api.jl")
elseif _is_patch(_GUROBI_VERSION, v"9.5")
    include("gen95/libgrb_common.jl")
    include("gen95/libgrb_api.jl")
elseif _is_patch(_GUROBI_VERSION, v"10.0")
    include("gen100/libgrb_common.jl")
    include("gen100/libgrb_api.jl")
elseif _is_patch(_GUROBI_VERSION, v"11.0")
    include("gen110/libgrb_common.jl")
    include("gen110/libgrb_api.jl")
else
    error("""
    You have installed version $_GUROBI_VERSION of Gurobi, which is not
    supported by Gurobi.jl. We require Gurobi version 9.0 or 9.1 or 9.5
    or 10.0 or 11.0.

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

for sym in filter(s -> startswith("$s", "GRB"), names(@__MODULE__, all = true))
    @eval export $sym
end

end
