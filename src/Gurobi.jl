# Copyright (c) 2015 Dahua Lin, Miles Lubin, Joey Huchette, Iain Dunning, and
#   contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

module Gurobi

# deps.jl file is always built via `Pkg.build`, even if we didn't find a local
# install and we want to use the artifact instead. This is so Gurobi.jl will be
# recompiled if we update the file. See issue #438 for more details.
const _DEPS_FILE = joinpath(dirname(@__FILE__), "..", "deps", "deps.jl")

if isfile(_DEPS_FILE)
    include(_DEPS_FILE)
else
    error(
        """
        Gurobi.jl is not installed correctly. Please run the following code and
        then restart Julia:
        ```
        import Pkg
        Pkg.build("Gurobi")
        ```
        """,
    )
end

if isdefined(@__MODULE__, :libgurobi)
    # deps.jl must define a local installation.
elseif Sys.islinux() || Sys.isapple() || Sys.iswindows()
    import Gurobi_jll: libgurobi
else
    error(
        "Unsupported platform: Use a manual installation by setting " *
        "`GUROBI_JL_USE_GUROBI_JLL` to false. See the README for details.",
    )
end

const _GUROBI_VERSION = if libgurobi == "__skipped_installation__"
    # The deps file is fake, with the intention to make Gurobi.jl loadable but
    # not usable.
    VersionNumber(10, 0, 0)
else
    let
        majorP, minorP, technicalP = Ref{Cint}(), Ref{Cint}(), Ref{Cint}()
        @ccall libgurobi.GRBversion(
            majorP::Ptr{Cint},
            minorP::Ptr{Cint},
            technicalP::Ptr{Cint},
        )::Cvoid
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
    include("gen100/libgrb_api.jl")
elseif _is_patch(_GUROBI_VERSION, v"11.0")
    include("gen110/libgrb_api.jl")
elseif _is_patch(_GUROBI_VERSION, v"12.0")
    include("gen120/libgrb_api.jl")
elseif _is_patch(_GUROBI_VERSION, v"13.0")
    include("gen130/libgrb_api.jl")
else
    error(
        """
        You have installed version $_GUROBI_VERSION of Gurobi, which is not
        supported by Gurobi.jl. We require Gurobi version 9.0 or 9.1 or 9.5
        or 10.0 or 11.0 or 12.0 or 13.0.

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
        """,
    )
end

include("MOI_wrapper/MOI_wrapper.jl")
include("MOI_wrapper/MOI_callbacks.jl")
include("MOI_wrapper/MOI_multi_objective.jl")
include("MOI_wrapper/MOI_indicator_constraint.jl")
include("MOI_wrapper/MOI_nonlinear.jl")

# Gurobi exports all `GRBXXX` symbols. If you don't want all of these symbols in
# your environment, then use `import Gurobi` instead of `using Gurobi`.

for sym in filter(s -> startswith("$s", "GRB"), names(@__MODULE__, all = true))
    @eval export $sym
end

import PrecompileTools

function _precompile_(env::Env)
    model = MOI.Utilities.CachingOptimizer(
        MOI.Utilities.UniversalFallback(MOI.Utilities.Model{Float64}()),
        MOI.instantiate(
            () -> Gurobi.Optimizer(env);
            with_bridge_type = Float64,
        ),
    )
    MOI.set(model, MOI.Silent(), true)
    x = MOI.add_variables(model, 3)
    MOI.supports(model, MOI.VariableName(), typeof(x[1]))
    MOI.set(model, MOI.VariableName(), x[1], "x1")
    MOI.set(model, MOI.VariablePrimalStart(), x[1], 0.0)
    MOI.add_constraint(model, x[1], MOI.ZeroOne())
    MOI.add_constraint(model, x[2], MOI.Integer())
    for F in (MOI.VariableIndex, MOI.ScalarAffineFunction{Float64})
        MOI.supports_constraint(model, F, MOI.GreaterThan{Float64})
        MOI.supports_constraint(model, F, MOI.LessThan{Float64})
        MOI.supports_constraint(model, F, MOI.EqualTo{Float64})
    end
    MOI.supports_constraint(model, MOI.VariableIndex, MOI.ZeroOne)
    MOI.supports_constraint(model, MOI.VariableIndex, MOI.Integer)
    MOI.add_constraint(model, x[1], MOI.GreaterThan(0.0))
    MOI.add_constraint(model, x[2], MOI.LessThan(0.0))
    MOI.add_constraint(model, x[3], MOI.EqualTo(0.0))
    MOI.add_constrained_variable(model, MOI.GreaterThan(0.0))
    MOI.add_constrained_variable(model, MOI.LessThan(0.0))
    MOI.add_constrained_variable(model, MOI.EqualTo(0.0))
    MOI.add_constrained_variable(model, MOI.Integer())
    MOI.add_constrained_variable(model, MOI.ZeroOne())
    set = (MOI.GreaterThan(0.0), MOI.LessThan(0.0))
    MOI.supports_add_constrained_variable(model, typeof(set))
    MOI.add_constrained_variable(model, set)
    f = 1.0 * x[1] + x[2] + x[3]
    c1 = MOI.add_constraint(model, f, MOI.GreaterThan(0.0))
    MOI.set(model, MOI.ConstraintName(), c1, "c1")
    MOI.supports(model, MOI.ConstraintName(), typeof(c1))
    MOI.add_constraint(model, f, MOI.LessThan(0.0))
    MOI.add_constraint(model, f, MOI.EqualTo(0.0))
    MOI.set(model, MOI.ObjectiveSense(), MOI.MAX_SENSE)
    MOI.supports(model, MOI.ObjectiveFunction{typeof(f)}())
    MOI.set(model, MOI.ObjectiveFunction{typeof(f)}(), f)
    MOI.optimize!(model)
    MOI.get(model, MOI.TerminationStatus())
    MOI.get(model, MOI.PrimalStatus())
    MOI.get(model, MOI.DualStatus())
    MOI.get(model, MOI.VariablePrimal(), x)
    return
end

PrecompileTools.@setup_workload begin
    PrecompileTools.@compile_workload begin
        try
            env = Env()
            _precompile_(env)
        catch
            nothing
        end
    end
end

end
