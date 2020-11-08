const _DEPRECATED_ERROR_MESSAGE = """
The C API of Gurobi.jl has been rewritten to expose the complete C API, and
all old functions have been removed. For more information, see the Discourse
announcement: https://discourse.julialang.org/t/ann-upcoming-breaking-changes-to-cplex-jl-and-gurobi-jl

Here is a brief summary of the changes.

* Constants have changed. For example `CB_MIPNODE` is now `GRB_CB_MIPNODE`
    to match the C API.
* Function names have changed. For example `free_env(env)` is now
    `GRBfreeenv(env)`.
* For users of `Gurobi.Optimizer()`, `model.inner` is now a pointer to the C
    model, instead of a `Gurobi.Model` object. However, conversion means that
    you should always pass `model` instead of `model.inner` to the low-level
    functions. For example:
    ```julia
    model = direct_model(Gurobi.Optimizer())
    grb_model = backend(model)  # grb_model is Gurobi.Optimizer
    # Old
    Gurobi.tune_model(grb_model.inner)
    # New
    GRBtunemodel(grb_model)
    ```
* Some functions have been removed entirely. For example:
    ```julia
    using JuMP, Gurobi
    model = direct_model(Gurobi.Optimizer())
    optimize!(model)
    grb_model = backend(model)
    stat = Gurobi.get_status_code(grb_model.inner)
    ```
    is now:
    ```julia
    using JuMP, Gurobi
    model = direct_model(Gurobi.Optimizer())
    optimize!(model)
    valueP = Ref{Cint}()
    grb_model = backend(model)
    ret = GRBgetintattr(grb_model, "Status", valueP)
    if ret != 0
        # Do something because the call failed
    end
    stat = valueP[]
    ```

The new API is more verbose, but the names and function arguments are now
identical to the C API, documentation for which is available at:
https://www.gurobi.com/documentation/9.0/refman/c_api_details.html

To revert to the old API, use:

    import Pkg
    Pkg.add(Pkg.PackageSpec(name = "Gurobi", version = v"0.8.1"))

Then restart Julia for the change to take effect.
"""
add_bvar!(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)
export add_bvar!

add_bvars!(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)
export add_bvars!

add_constr!(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)
export add_constr!

add_constrs!(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)
export add_constrs!

add_constrs_t!(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)
export add_constrs_t!

add_cvar!(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)
export add_cvar!

add_cvars!(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)
export add_cvars!

add_diag_qpterms!(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)

add_ivar!(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)
export add_ivar!

add_ivars!(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)
export add_ivars!

add_qconstr!(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)
export add_qconstr!

add_qpterms!(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)
export add_qpterms!

add_rangeconstr!(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)
export add_rangeconstr!

add_rangeconstrs!(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)
export add_rangeconstrs!

add_rangeconstrs_t!(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)
export add_rangeconstrs_t!

add_sos!(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)
export add_sos!

add_var!(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)
export add_var!

add_vars!(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)
export add_vars!

c_api_setobjectiven(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)

cbcut(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)
export cbcut

cbget(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)

cbget_barrier_compl(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)
export cbget_barrier_compl

cbget_barrier_dualinf(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)
export cbget_barrier_dualinf

cbget_barrier_dualobj(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)
export cbget_barrier_dualobj

cbget_barrier_itrcnt(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)
export cbget_barrier_itrcnt

cbget_barrier_priminf(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)
export cbget_barrier_priminf

cbget_barrier_primobj(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)
export cbget_barrier_primobj

cbget_mip_cutcnt(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)
export cbget_mip_cutcnt

cbget_mip_itrcnt(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)
export cbget_mip_itrcnt

cbget_mip_nodcnt(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)
export cbget_mip_nodcnt

cbget_mip_nodlft(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)
export cbget_mip_nodlft

cbget_mip_objbnd(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)
export cbget_mip_objbnd

cbget_mip_objbst(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)
export cbget_mip_objbst

cbget_mip_solcnt(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)
export cbget_mip_solcnt

cbget_mipnode_nodcnt(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)
export cbget_mipnode_nodcnt

cbget_mipnode_objbnd(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)
export cbget_mipnode_objbnd

cbget_mipnode_objbst(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)
export cbget_mipnode_objbst

cbget_mipnode_rel(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)
export cbget_mipnode_rel

cbget_mipnode_solcnt(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)
export cbget_mipnode_solcnt

cbget_mipnode_status(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)
export cbget_mipnode_status

cbget_mipsol_nodcnt(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)
export cbget_mipsol_nodcnt

cbget_mipsol_obj(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)
export cbget_mipsol_obj

cbget_mipsol_objbnd(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)
export cbget_mipsol_objbnd

cbget_mipsol_objbst(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)
export cbget_mipsol_objbst

cbget_mipsol_rel(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)

cbget_mipsol_sol(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)
export cbget_mipsol_sol

cbget_mipsol_solcnt(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)
export cbget_mipsol_solcnt

cbget_pre_bndchg(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)
export cbget_pre_bndchg

cbget_pre_coechg(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)
export cbget_pre_coechg

cbget_pre_coldel(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)
export cbget_pre_coldel

cbget_pre_rowdel(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)
export cbget_pre_rowdel

cbget_pre_senchg(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)
export cbget_pre_senchg

cbget_runtime(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)
export cbget_runtime

cbget_spx_dualinf(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)
export cbget_spx_dualinf

cbget_spx_ispert(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)
export cbget_spx_ispert

cbget_spx_itrcnt(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)
export cbget_spx_itrcnt

cbget_spx_objval(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)
export cbget_spx_objval

cbget_spx_priminf(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)
export cbget_spx_priminf

cblazy(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)
export cblazy

cbsolution(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)

cchar(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)

check_moi_callback_validity(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)

checkvalue(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)

chg_coeffs!(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)
export chg_coeffs!

computeIIS(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)
export computeIIS

compute_conflict(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)

copy(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)
# export copy # Don't export this copy.

cvec(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)

cvecx(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)

default_moi_callback(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)

del_constrs!(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)
export del_constrs!

del_sos!(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)

del_vars!(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)
export del_vars!

delq!(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)

delqconstrs!(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)

eval(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)

fixed_model(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)
export fixed_model

free_env(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)
export free_env

free_model(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)

fvec(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)

fvecx(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)

get_barrier_iter_count(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)

get_basis(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)

get_charattrarray(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)

get_charattrarray!(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)

get_charattrelement(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)

get_charattrlist(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)

get_charattrlist!(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)

get_constrmatrix(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)
export get_constrmatrix

get_constrs(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)

get_dbl_param(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)

get_dblattr(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)

get_dblattrarray(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)

get_dblattrarray!(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)

get_dblattrelement(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)

get_dblattrlist(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)

get_dblattrlist!(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)

get_error_msg(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)

get_int_param(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)

get_intattr(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)

get_intattrarray(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)

get_intattrarray!(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)

get_intattrelement(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)

get_intattrlist(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)

get_intattrlist!(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)

get_iter_count(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)

get_multiobj_c(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)

get_multiobj_n(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)

get_multiobj_priority(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)

get_multiobj_weight(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)

get_node_count(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)

get_objbound(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)

get_objval(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)
export get_objval

get_optiminfo(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)
export get_optiminfo

get_runtime(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)

get_sol_count(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)

get_solution(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)
export get_solution

get_sos(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)

get_sos_matrix(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)

get_status(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)
export get_status

get_status_code(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)

get_str_param(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)

get_strattr(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)

get_strattrarray(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)

get_strattrarray!(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)

get_strattrelement(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)

get_tune_result!(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)
export get_tune_result!

getcoeff(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)

getcoeff!(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)

getlibversion(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)

getparam(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)
export getparam

getq(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)

getqconstr(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)

gurobi_callback_wrapper(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)

gurobi_model(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)
export gurobi_model

include(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)

is_mip(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)
export is_mip

is_qcp(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)
export is_qcp

is_qp(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)
export is_qp

is_valid(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)

ivec(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)

loadbasis(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)

lowerbounds(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)
export lowerbounds

mastercallback(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)

model_name(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)
export model_name

model_sense(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)
export model_sense

model_type(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)
export model_type

num_binvars(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)

num_cnzs(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)
export num_cnzs

num_constrs(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)
export num_constrs

num_intvars(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)

num_qcnzs(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)
export num_qcnzs

num_qconstrs(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)
export num_qconstrs

num_qnzs(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)
export num_qnzs

num_sos(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)
export num_sos

num_vars(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)
export num_vars

objcoeffs(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)
export objcoeffs

optimize(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)
export optimize

presolve_model(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)
export presolve_model

read(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)
# export read # Don't export this read

read_model(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)
export read_model

reset_model!(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)
export reset_model!

set_callback_func!(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)
export set_callback_func!

set_charattrarray!(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)

set_charattrelement!(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)

set_charattrlist!(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)

set_dbl_param!(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)

set_dblattr!(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)

set_dblattrarray!(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)

set_dblattrelement!(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)

set_dblattrlist!(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)

set_int_param!(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)

set_intattr!(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)

set_intattrarray!(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)

set_intattrelement!(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)

set_intattrlist!(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)

set_multiobj!(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)

set_multiobj_c!(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)

set_multiobj_n!(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)

set_multiobj_priority!(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)

set_multiobj_weight!(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)

set_objcoeffs!(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)
export set_objcoeffs!

set_sense!(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)
export set_sense!

set_str_param!(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)

set_strattr!(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)

set_strattrelement!(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)

setmathprogcallback!(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)

setparam!(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)
export setparam!

setparams!(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)
export setparams!

sparse_transpose(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)

terminate(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)

tune_model(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)
export tune_model

update_model!(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)
export update_model!

updatemodel!(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)

upperbounds(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)
export upperbounds

write_model(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)
export write_model

GurobiSolver(args...; kwargs...) = error(_DEPRECATED_ERROR_MESSAGE)
export GurobiSolver

