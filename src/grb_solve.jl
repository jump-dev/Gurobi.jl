# Optimization and solution query

function optimize(model::Model)
    @assert model.ptr_model != C_NULL
    ret = @grb_ccall(optimize, Cint, (Ptr{Void},), model)
    if ret != 0
        throw(GurobiError(model.env, ret))
    end
    nothing
end

function computeIIS(model::Model)
    @assert model.ptr_model != C_NULL
    ret = @grb_ccall(computeIIS, Cint, (Ptr{Void},), model)
    if ret != 0
        throw(GurobiError(model.env, ret))
    end
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
    :suboptimal
]

get_status_code(model::Model) = get_intattr(model, "Status")
get_status(model::Model) = status_symbols[get_status_code(model)]::Symbol

@grb_dbl_attr get_runtime      "Runtime"
@grb_dbl_attr get_objval       "ObjVal"
@grb_dbl_attr get_objbound     "ObjBound"

@grb_int_attr get_sol_count    "SolCount"
@grb_int_attr get_barrier_iter_count "BarIterCount"
@grb_dbl_attr get_node_count   "NodeCount"

get_iter_count(model::Model) = convert(Int, get_dblattr(model, "IterCount"))
get_node_count(model::Model) = convert(Int, get_dblattr(model, "NodeCount"))


type OptimInfo
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
    -1 => :Nonbasic)

function get_basis(model::Model)
    cval = Array(Cint, num_vars(model))
    cbasis = Array(Symbol, num_vars(model))
    get_intattrarray!(cval, model, "VBasis", 1)
    for it in 1:length(cval)
        cbasis[it] = varmap[cval[it]]
    end
    
    rval = Array(Cint, num_constrs(model))
    rbasis = Array(Symbol, num_constrs(model))
    get_intattrarray!(rval, model, "CBasis", 1)
    rsense = Array(Cchar, num_constrs(model))
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
