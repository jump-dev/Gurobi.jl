module Gurobi

const _DEPS_FILE = joinpath(dirname(@__FILE__),"..","deps","deps.jl")
if isfile(_DEPS_FILE)
    include(_DEPS_FILE)
else
    error("""
        Gurobi not properly installed. Please run Pkg.build(\"Gurobi\"). For
        more information go to https://github.com/JuliaOpt/Gurobi.jl
    """)
end

import SparseArrays, LinearAlgebra

### imports

import Base.show, Base.copy, Base.read

### exports
export

# grb_env
free_env,

# grb_params
getparam, setparam!, setparams!,

# grb_model
set_sense!, update_model!, reset_model!, get_tune_result!,
read_model, write_model, tune_model, presolve_model, fixed_model,
copy, read,

# grb_attrs
model_name, model_sense, model_type,
num_vars, num_constrs, num_sos, num_qconstrs,
num_cnzs, num_qnzs, num_qcnzs,
is_qp, is_qcp, is_mip,

lowerbounds, upperbounds, objcoeffs, set_objcoeffs!,

# grb_vars
GRB_CONTINUOUS, GRB_BINARY, GRB_INTEGER,
add_var!, add_vars!, add_cvar!, add_cvars!,
add_bvar!, add_bvars!, add_ivar!, add_ivars!,
del_vars!,

# grb_constrs
add_constr!, add_constrs!, add_constrs_t!,
add_rangeconstr!, add_rangeconstrs!, add_rangeconstrs_t!,
get_constrmatrix, add_sos!, del_constrs!, chg_coeffs!,

# grb_quad
add_qpterms!, add_qconstr!,

# higher level
gurobi_model,

# grb_solve
optimize, computeIIS, get_solution,
get_status, OptimInfo, get_optiminfo, get_objval

### include source files

include("grb_common.jl")
include("grb_env.jl")

include("grb_model.jl")
include("grb_params.jl")
include("grb_vars.jl")
include("grb_attrs.jl")
include("grb_constrs.jl")
include("grb_quad.jl")
include("grb_highlevel.jl")

include("grb_solve.jl")
include("grb_callbacks.jl")

include("MPB_wrapper.jl")
include("MOI_wrapper.jl")

end
