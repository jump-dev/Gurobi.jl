# Optimization and solution query

function optimize(model::Model)
    @assert model.ptr_model != C_NULL
    ret = @grb_ccall(optimize, Cint, (Ptr{Cvoid},), model)
    if ret != 0
        throw(GurobiError(model.env, ret))
    end
    nothing
end

function computeIIS(model::Model)
    @assert model.ptr_model != C_NULL
    ret = @grb_ccall(computeIIS, Cint, (Ptr{Cvoid},), model)
    if ret != 0
        throw(GurobiError(model.env, ret))
    end
    model.conflict = ret
    nothing
end


#################################################
#
#  solution status and optimization info
#
#################################################

const GRB_LOADED          = 1
const GRB_OPTIMAL         = 2
const GRB_INFEASIBLE      = 3
const GRB_INF_OR_UNBD     = 4
const GRB_BOUNDED         = 5
const GRB_CUTOFF          = 6
const GRB_ITERATION_LIMIT = 7
const GRB_NODE_LIMIT      = 8
const GRB_TIME_LIMIT      = 9
const GRB_SOLUTION_LIMIT  = 10
const GRB_INTERRUPTED     = 11
const GRB_NUMERIC         = 12
const GRB_SUBOPTIMAL      = 13
const GRB_INPROGRESS      = 14
const GRB_USER_OBJ_LIMIT  = 15

const status_symbols = [
    :loaded,
    :optimal,
    :infeasible,
    :inf_or_unbd,
    :unbounded,
    :cutoff,
    :iteration_limit,
    :node_limit,
    :time_limit,
    :solution_limit,
    :interrupted,
    :numeric,
    :suboptimal,
    :inprogress,
    :user_obj_limit
]

get_status_code(model::Model) = get_intattr(model, "Status")
get_status(model::Model) = status_symbols[get_status_code(model)]::Symbol

@grb_dbl_attr get_runtime      "Runtime"
@grb_dbl_attr get_objval       "ObjVal"
@grb_dbl_attr get_objbound     "ObjBound"

@grb_int_attr get_sol_count    "SolCount"
@grb_int_attr get_barrier_iter_count "BarIterCount"
#@grb_dbl_attr get_node_count   "NodeCount"

get_iter_count(model::Model) = convert(Int, get_dblattr(model, "IterCount"))
get_node_count(model::Model) = convert(Int, get_dblattr(model, "NodeCount"))


mutable struct OptimInfo
    status::Symbol
    runtime::Float64

    sol_count::Int
    iter_count::Int
    barrier_iter_count::Int
    node_count::Int
end

function get_optiminfo(model::Model)
    OptimInfo(
        get_status(model),
        get_runtime(model),

        get_sol_count(model),
        get_iter_count(model),
        get_barrier_iter_count(model),
        get_node_count(model)
    )
end

function show(io::IO, s::OptimInfo)
    println(io, "Gurobi Optimization Info")
    println(io, "    status   = $(s.status)")
    println(io, "    runtime  = $(s.runtime)")
    println(io, "    # solutions = $(s.sol_count)")
    println(io, "    # iters     = $(s.iter_count)")
    println(io, "    # bar iters = $(s.barrier_iter_count)")
    println(io, "    # nodes     = $(s.node_count)")
end

#################################################
#
#  solution query
#
#################################################

get_solution(model::Model) = get_dblattrarray(model, "X", 1, num_vars(model))

const varmap = Dict(
    -3 => :Superbasic,
    -2 => :NonbasicAtUpper,
    -1 => :NonbasicAtLower,
     0 => :Basic
)

const conmap = Dict(
     0 => :Basic,
    -1 => :Nonbasic
)

const basicmap_rev = Dict(
    :Superbasic => Cint(-3),
    :NonbasicAtUpper => Cint(-2),
    :NonbasicAtLower => Cint(-1),
    :Basic => Cint(0),
    :Nonbasic => Cint(-1)
)

function get_basis(model::Model)
    cval = Array{Cint}(undef, num_vars(model))
    cbasis = Array{Symbol}(undef, num_vars(model))
    get_intattrarray!(cval, model, "VBasis", 1)
    for it in 1:length(cval)
        cbasis[it] = varmap[cval[it]]
    end

    rval = Array{Cint}(undef, num_constrs(model))
    rbasis = Array{Symbol}(undef, num_constrs(model))
    get_intattrarray!(rval, model, "CBasis", 1)
    rsense = Array{Cchar}(undef, num_constrs(model))
    get_charattrarray!(rsense, model, "Sense", 1)
    for it in 1:length(rval)
        rbasis[it] = conmap[rval[it]]
        if rbasis[it] == :Nonbasic
            if rsense[it] == convert(Cchar,'<')
                rbasis[it] = :NonbasicAtUpper
            else
                rbasis[it] = :NonbasicAtLower
            end
        end
    end
    return cbasis, rbasis
end

"""
    loadbasis(model::Model, x::Vector)

Load basis to a problem in form of primal solution.

    loadbasis(model::Model, rval::Vector{Symbol}, cval::Vector{Symbol})

Load basis to a problem in terms of basicness description Variables columns (in
`cval`) can be: `:Basic`, `:NonbasicAtLower`, `:NonbasicAtUpper` and
`:Superbasic` Constraints rows (in `rval`) can be: `:Basic`, `:Nonbasic`
"""
function loadbasis(model::Model, x::Vector)

    ncols = num_vars(model)
    nrows = num_constrs(model)

    length(x) != ncols && error("solution candidate size is different from the number of columns")

    cvals = Array{Cint}(undef, ncols)
    rvals = Array{Cint}(undef, nrows)

    # obtain situation of columns

    lb = get_dblattrarray(model, "LB", 1, ncols)
    ub = get_dblattrarray(model, "UB", 1, ncols)

    for i in 1:ncols
        if isapprox(x[i],lb[i])
            cvals[i] = basicmap_rev[:NonbasicAtLower]
        elseif isapprox(x[i],ub[i])
            cvals[i] = basicmap_rev[:NonbasicAtUpper]
        else
            cvals[i] = basicmap_rev[:Basic]
        end
    end

    # obtain situation of rows where: y = Ax

    A = get_constrmatrix(model) #A is sparse

    y = A*x

    senses = get_charattrarray(model, "Sense", 1, nrows)
    rhs    = get_dblattrarray(model, "RHS", 1, nrows)

    for j in 1:nrows
        if senses[j] == '=' && isapprox(y[j], rhs[j])
            rvals[j] = basicmap_rev[:Nonbasic]#AtLower
        elseif senses[j] == '>' && isapprox(y[j], rhs[j])
            rvals[j] = basicmap_rev[:Nonbasic]#AtLower
        elseif senses[j] == '<' && isapprox(y[j], rhs[j])
            rvals[j] = basicmap_rev[:Nonbasic]#AtUpper
        else
            rvals[j] = basicmap_rev[:Basic]
        end
    end

    loadbasis(model, rvals, cvals)

    return nothing
end
function loadbasis(model::Model, rval::Vector{Symbol}, cval::Vector{Symbol})

    nrval = map(x->conmap[x], rval)
    ncval = map(x->varmap[x], cval)

    loadbasis(model, nrval, ncval)

    return nothing
end
function loadbasis(model::Model, rval::Vector{Cint}, cval::Vector{Cint})

    ncols = num_vars(model)
    nrows = num_constrs(model)

    set_intattrarray!(model, "VBasis", 1, num_vars(model), cval)
    set_intattrarray!(model, "CBasis", 1, num_constrs(model), rval) # r = row; c = constr

    return nothing
end
