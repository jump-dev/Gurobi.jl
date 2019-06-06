# GurobiMathProgInterface
# Standardized MILP interface

export GurobiSolver

import MathProgBase
const MPB = MathProgBase

mutable struct GurobiMathProgModel <: MPB.AbstractLinearQuadraticModel
    inner::Model
    last_op_type::Symbol  # To support arbitrary order of addVar/addCon
                          # Two possibilities :Var :Con
    changed_constr_bounds::Bool # have we updated the bounds below?
    obj::Vector{Float64} # objective vector
    lb::Vector{Float64} # persistent bounds on constraints to maintain
    ub::Vector{Float64} # abstraction for lb ≤ Ax ≤ ub
    lazycb
    cutcb
    heuristiccb
    infocb
    options
end

function GurobiMathProgModel(env=nothing;options...)
   finalize_env = (env == nothing)
   if env == nothing
       env = Env()
   end
   m = GurobiMathProgModel(Model(env,""; finalize_env=finalize_env), :Con, false, Float64[], Float64[], Float64[], nothing, nothing, nothing, nothing, options)
   setparams!(m)
   return m
end

function Base.copy(m::GurobiMathProgModel)

    m.lazycb == nothing || @warn("Callbacks can't be copied, lazy callback ignored")
    m.cutcb == nothing || @warn("Callbacks can't be copied, cut callback ignored")
    m.heuristiccb == nothing || @warn("Callbacks can't be copied, heuristic callback ignored")
    m.infocb == nothing || @warn("Callbacks can't be copied, info callback ignored")

    return GurobiMathProgModel(copy(m.inner),
                               m.last_op_type,
                               m.changed_constr_bounds,
                               copy(m.obj),
                               copy(m.lb),
                               copy(m.ub),
                               nothing,
                               nothing,
                               nothing,
                               nothing,
                               deepcopy(m.options))
end

function setparams!(m::GurobiMathProgModel)
    # Helper to set the parameters on the model's copy of env rather than
    # modifying the global env (ref: http://www.gurobi.com/support/faqs#P)

    # Set `InfUnbdInfo` to 1 by default so infeasibility rays available
    setparam!(m.inner, "InfUnbdInfo", 1)
    for (name,value) in m.options
        setparam!(m.inner, string(name), value)
    end
end


mutable struct GurobiSolver <: MPB.AbstractMathProgSolver
    env
    options
end
GurobiSolver(env=nothing; kwargs...) = GurobiSolver(env, kwargs)
MPB.LinearQuadraticModel(s::GurobiSolver) = GurobiMathProgModel(s.env; s.options...)
MPB.ConicModel(s::GurobiSolver) = MPB.LPQPtoConicBridge(MPB.LinearQuadraticModel(s))

MPB.supportedcones(::GurobiSolver) = [:Free,:Zero,:NonNeg,:NonPos,:SOC,:SOCRotated]

function MPB.setparameters!(s::GurobiSolver; mpboptions...)
    opts = collect(Any,s.options)
    for (optname, optval) in mpboptions
        if optname == :TimeLimit
            push!(opts, (:TimeLimit,optval))
        elseif optname == :Silent
            if optval == true
                push!(opts, (:OutputFlag,0))
            end
        else
            error("Unrecognized parameter $optname")
        end
    end
    s.options = opts
    return
end

function MPB.setparameters!(m::GurobiMathProgModel; mpboptions...)
    for (optname, optval) in mpboptions
        if optname == :TimeLimit
            setparam!(m.inner, "TimeLimit", optval)
        elseif optname == :Silent
            if optval == true
                setparam!(m.inner,"OutputFlag",0)
            end
        else
            error("Unrecognized parameter $optname")
        end
    end
end

function MPB.loadproblem!(m::GurobiMathProgModel, filename::AbstractString)
    read_model(m.inner, filename)
    m.obj = MPB.getobj(m)
    if checkvalue(m.obj, GRB_INFINITY)
        _objwarning(m.obj)
    end
    _truncateobj!(m.obj)
    m.lb = MPB.getconstrLB(m)
    m.ub = MPB.getconstrUB(m)
end

function MPB.loadproblem!(m::GurobiMathProgModel, A::LinearAlgebra.Adjoint{T, Array{T, 2}},
                        collb, colub, obj, rowlb, rowub, sense) where T
    MPB.loadproblem!(m, collect(A), collb, colub, obj, rowlb, rowub, sense)
end

function MPB.loadproblem!(m::GurobiMathProgModel, A, collb, colub, obj, rowlb, rowub, sense)
  # throw away old model but keep env and finalize_env
  env = m.inner.env
  finalize_env = m.inner.finalize_env
  m.inner.finalize_env = false
  free_model(m.inner)
  m.inner = Model(env, "", finalize_env=finalize_env)
  # Re-set options on new model's copy of env
  setparams!(m)

  add_cvars!(m.inner, float(copy(obj)), float(collb), float(colub))
  update_model!(m.inner)

  neginf = typemin(eltype(rowlb))
  posinf = typemax(eltype(rowub))

  # check if we have any range constraints
  # to properly support these, we will need to keep track of the
  # slack variables automatically added by gurobi.
  rangeconstrs = sum((rowlb .!= rowub) .& (rowlb .> neginf) .& (rowub .< posinf))
  if rangeconstrs > 0
      @warn("Julia Gurobi interface doesn't properly support range " *
            "(two-sided) constraints. See Gurobi.jl issue #14")
      add_rangeconstrs!(m.inner, float(A), float(rowlb), float(rowub))
      # A work around for the additional slack variables introduced and the
      # reformulation bug warning
      append!(obj, zeros(rangeconstrs))
  else
      b = Array{Float64}(undef, length(rowlb))
      senses = Array{Cchar}(undef, length(rowlb))
      for i in 1:length(rowlb)
          if rowlb[i] == rowub[i]
              senses[i] = '='
              b[i] = rowlb[i]
          elseif rowlb[i] > neginf
              senses[i] = '>'
              b[i] = rowlb[i]
          else
              @assert rowub[i] < posinf
              senses[i] = '<'
              b[i] = rowub[i]
          end
      end
      add_constrs!(m.inner, float(A), senses, b)
  end

  m.lb, m.ub = rowlb, rowub
  m.obj = copy(obj)
  _truncateobj!(m.obj)
  update_model!(m.inner)
  MPB.setsense!(m,sense)
end

function _truncateobj(v::Real)
    # Gurobi truncates objective coefficients to +/-GRB_INFINITY
    # We need to do this to the m.obj vector for consistency
    #   see #94
    if v > GRB_INFINITY
        return GRB_INFINITY
    elseif v < -GRB_INFINITY
        return -GRB_INFINITY
    end
    return v
end
_truncateobj!(obj::Vector) = map!(_truncateobj, obj, obj)

MPB.writeproblem(m::GurobiMathProgModel, filename::AbstractString) = write_model(m.inner, filename)

MPB.getvarLB(m::GurobiMathProgModel)     = get_dblattrarray( m.inner, "LB", 1, num_vars(m.inner))
function MPB.setvarLB!(m::GurobiMathProgModel, l)
    if checkvalue(l, GRB_BOUNDMAX)
        _boundwarning(l, getvarUB(m))
    end
    set_dblattrarray!(m.inner, "LB", 1, num_vars(m.inner), l)
end

MPB.getvarUB(m::GurobiMathProgModel)     = get_dblattrarray( m.inner, "UB", 1, num_vars(m.inner))
function MPB.setvarUB!(m::GurobiMathProgModel, u)
    if checkvalue(u, GRB_BOUNDMAX)
        _boundwarning(getvarLB(m), u)
    end
    set_dblattrarray!(m.inner, "UB", 1, num_vars(m.inner), u)
end

function MPB.getconstrLB(m::GurobiMathProgModel)
    sense = get_charattrarray(m.inner, "Sense", 1, num_constrs(m.inner))
    ret   = get_dblattrarray(m.inner, "RHS", 1, num_constrs(m.inner))
    for i = 1:num_constrs(m.inner)
        if sense[i] == '>' || sense[i] == '='
            # Do nothing
        else
            # LEQ constraint, so LB is -Inf
            ret[i] = -Inf
        end
     end
     return ret
end

function MPB.getconstrUB(m::GurobiMathProgModel)
    sense = get_charattrarray(m.inner, "Sense", 1, num_constrs(m.inner))
    ret   = get_dblattrarray(m.inner, "RHS", 1, num_constrs(m.inner))
    for i = 1:num_constrs(m.inner)
        if sense[i] == '<' || sense[i] == '='
            # Do nothing
        else
            # GEQ constraint, so UB is +Inf
            ret[i] = +Inf
        end
    end
    return ret
end

# setconstrLB!(m::GurobiMathProgModel, lb) = (m.changed_constr_bounds = true; m.last_op_type = :Con; m.lb = copy(lb))
# setconstrUB!(m::GurobiMathProgModel, ub) = (m.changed_constr_bounds = true; m.last_op_type = :Con; m.ub = copy(ub))
MPB.setconstrLB!(m::GurobiMathProgModel, lb) = (m.changed_constr_bounds = true; m.lb = copy(lb))
MPB.setconstrUB!(m::GurobiMathProgModel, ub) = (m.changed_constr_bounds = true; m.ub = copy(ub))

MPB.getobj(m::GurobiMathProgModel, i::Int=1) = get_dblattrarray( m.inner, "Obj", i, num_vars(m.inner))
function MPB.setobj!(m::GurobiMathProgModel, c, i::Int=1)
    if checkvalue(c, GRB_INFINITY)
        _objwarning(c)
    end
    m.obj = copy(c)
    _truncateobj!(m.obj)
    set_dblattrarray!(m.inner, "Obj", i, num_vars(m.inner), c)
end

set_multiobj_n!(m::Gurobi.Model, n::Int) = Gurobi.set_intattr!(m, "NumObj", n)
get_multiobj_n(m::Gurobi.Model)          = Gurobi.get_intattr(m, "NumObj")

function set_multiobj_c!(m::Model, i::Int, c::AbstractVector{Float64})
    _chklen(c, num_vars(m))
    i0 = get_int_param(m, "ObjNumber")
    set_int_param!(m, "ObjNumber", i-1)
    set_dblattrarray!(m, "ObjN", 1, num_vars(m), c)
    set_int_param!(m, "ObjNumber", i0)
end

function get_multiobj_c(m::Model, i::Int)
    i0 = get_int_param(m, "ObjNumber")
    set_int_param!(m, "ObjNumber", i-1)
    c = get_dblattrarray(m, "ObjN", 1, num_vars(m))
    set_int_param!(m, "ObjNumber", i0)
    return c
end

function set_multiobj_priority!(m::Model, i::Int, priority::Int)
    i0 = get_int_param(m, "ObjNumber")
    set_int_param!(m, "ObjNumber", i-1)
    set_intattr!(m, "ObjNPriority", priority)
    set_int_param!(m, "ObjNumber", i0)
end

function get_multiobj_priority(m::Model, i::Int)::Int
    i0 = get_int_param(m, "ObjNumber")
    set_int_param!(m, "ObjNumber", i-1)
    p = get_intattr(m, "ObjNPriority")
    set_int_param!(m, "ObjNumber", i0)
    return p
end

function set_multiobj_weight!(m::Model, i::Int, w::Float64)
    i0 = get_int_param(m, "ObjNumber")
    set_int_param!(m, "ObjNumber", i-1)
    set_dblattr!(m, "ObjNWeight", w)
    set_int_param!(m, "ObjNumber", i0)
end

function get_multiobj_weight(m::Model, i::Int)::Float64
    i0 = get_int_param(m, "ObjNumber")
    set_int_param!(m, "ObjNumber", i-1)
    w::Float64 = get_dblattr(m, "ObjNWeight")
    set_int_param!(m, "ObjNumber", i0)
    return w
end

function set_multiobj!(m::Model, i::Int, c::AbstractVector{Float64}, p::Int, w::Float64)
    nobj = get_multiobj_n(m)
    if !(1 ≤ i ≤ nobj)
        error("Tried to set objective $i of $nobj. You must use the set_multiobj_n to specify the correct number of objective functions for your problem.")
    end
    set_multiobj_c!(m, i, c)
    set_multiobj_priority!(m, i, p)
    set_multiobj_weight!(m, i, w)
end

function MPB.addvar!(m::GurobiMathProgModel, constridx, constrcoef, l, u, objcoef)
    if m.last_op_type == :Con
        updatemodel!(m)
        m.last_op_type = :Var
    end
    push!(m.obj, _truncateobj(objcoef))
    add_var!(m.inner, length(constridx), constridx, float(constrcoef), float(objcoef), float(l), float(u), GRB_CONTINUOUS)
end
MPB.addvar!(m::GurobiMathProgModel, l, u, objcoef) = MPB.addvar!(m, Int[], Float64[], l, u, objcoef)

function MPB.delvars!(m::GurobiMathProgModel, idx)
    deleteat!(m.obj, idx)
    del_vars!(m.inner, idx)
end

function MPB.addconstr!(m::GurobiMathProgModel, varidx, coef, lb, ub)
    if m.last_op_type == :Var
        updatemodel!(m)
        m.last_op_type = :Con
    end
    if lb == -Inf
        # <= constraint
        add_constr!(m.inner, varidx, coef, '<', ub)
    elseif ub == +Inf
        # >= constraint
        add_constr!(m.inner, varidx, coef, '>', lb)
    elseif lb == ub
        # == constraint
        add_constr!(m.inner, varidx, coef, '=', lb)
    else
        # Range constraint
        error("Adding range constraints not supported yet.")
    end
end
MPB.delconstrs!(m::GurobiMathProgModel, idx) = del_constrs!(m.inner, idx)

MPB.changecoeffs!(m::GurobiMathProgModel, cidx, vidx, val) = chg_coeffs!(m.inner, cidx, vidx, val)

function updatemodel!(m::GurobiMathProgModel)
    update_model!(m.inner)
    if Gurobi.version < v"7.0" && m.obj != MPB.getobj(m)
        @warn("""
            You have encountered a known bug in Gurobi. Any information you query from the model may be incorrect.
            This bug has existed since the first version of Gurobi but is fixed in Gurobi v7.0.

            For more information go to https://github.com/JuliaOpt/Gurobi.jl/issues/60
            Please leave a comment stating that you encountered this bug! We would like to know how prevalent it is.
        """)
    end
    if m.changed_constr_bounds
    # update lower/upper bounds on linear constraints (if they're consistent...)
        sense = get_charattrarray(m.inner, "Sense", 1, num_constrs(m.inner))
        rhs   = get_dblattrarray( m.inner, "RHS",   1, num_constrs(m.inner))
        lb, ub = m.lb, m.ub
        @assert (n_rows = num_constrs(m.inner)) == length(lb) == length(ub) ==
                                                   length(sense) == length(rhs)
        for i = 1:n_rows
            if -Inf < lb[i] == ub[i] < Inf
                sense[i] = '='
                rhs[i]   = lb[i]
            elseif (-Inf < lb[i] && ub[i] < Inf)
                error("Gurobi.jl does not currently support range constraints")
            elseif (lb[i] == -Inf && ub[i] < Inf)
                sense[i] = '<'
                rhs[i]   = ub[i]
            elseif (lb[i] > -Inf && ub[i] == Inf)
                sense[i] = '>'
                rhs[i]   = lb[i]
            else
                error("Internal error.")
            end
        end
        set_charattrarray!(m.inner, "Sense", 1, num_constrs(m.inner), sense)
        set_dblattrarray!(m.inner, "RHS", 1, num_constrs(m.inner), rhs)
        m.changed_constr_bounds = false
        update_model!(m.inner)
    end
end

MPB.getconstrmatrix(m::GurobiMathProgModel) = get_constrmatrix(m.inner)

function MPB.setsense!(m::GurobiMathProgModel, sense)
  if sense == :Min
    set_sense!(m.inner, :minimize)
  elseif sense == :Max
    set_sense!(m.inner, :maximize)
  else
    error("Unrecognized objective sense $sense")
  end
end
function MPB.getsense(m::GurobiMathProgModel)
  v = get_intattr(m.inner, "ModelSense")
  if v == -1
    return :Max
  else
    return :Min
  end
end

MPB.numvar(m::GurobiMathProgModel)    = (updatemodel!(m); num_vars(m.inner))
MPB.numconstr(m::GurobiMathProgModel) = num_constrs(m.inner) + num_qconstrs(m.inner)
MPB.numlinconstr(m::GurobiMathProgModel) = num_constrs(m.inner)
MPB.numquadconstr(m::GurobiMathProgModel) = num_qconstrs(m.inner)

function MPB.optimize!(m::GurobiMathProgModel)
    # set callbacks if present
    if m.lazycb != nothing || m.cutcb != nothing || m.heuristiccb != nothing || m.infocb != nothing
        updatemodel!(m)
        if !is_mip(m.inner)
            @warn("Gurobi ignores branch-and-bound callbacks when no discrete elements are present in the model.")
        end
        setmathprogcallback!(m)
    end
    if m.lazycb != nothing
      setparam!(m.inner, "LazyConstraints", 1)
    end
    updatemodel!(m)
    optimize(m.inner)
end

function MPB.status(m::GurobiMathProgModel)
  s = get_status(m.inner)
  if s == :optimal
    return :Optimal
  elseif s == :infeasible
    return :Infeasible
  elseif s == :unbounded
    return :Unbounded
  elseif s == :inf_or_unbd
    @warn("Gurobi reported infeasible or unbounded. Set InfUnbdInfo=1 for more specific status.")
    return :InfeasibleOrUnbounded
  elseif s == :iteration_limit || s == :node_limit || s == :time_limit || s == :solution_limit
    return :UserLimit
  elseif s == :numeric
    return :Numeric
  elseif s == :suboptimal
    return :Suboptimal # not very useful status
  elseif s == :interrupted # ended by user?
    return :UserLimit
  elseif s == :inprogress
    return :InProgress
  elseif s == :user_obj_limit
    return :UserObjLimit
  elseif s == :cutoff
    return :Cutoff
  else
    error("Unrecognized solution status: $s")
  end
end

MPB.getobjval(m::GurobiMathProgModel)   = get_objval(m.inner)
MPB.getobjbound(m::GurobiMathProgModel) = get_objbound(m.inner)
MPB.getsolution(m::GurobiMathProgModel) = get_solution(m.inner)

function MPB.getconstrsolution(m::GurobiMathProgModel)
    sense = get_charattrarray(m.inner, "Sense", 1, num_constrs(m.inner))
    rhs   = get_dblattrarray( m.inner, "RHS",   1, num_constrs(m.inner))
    slack = get_dblattrarray( m.inner, "Slack", 1, num_constrs(m.inner))
    ret   = zeros(num_constrs(m.inner))
    for i = 1:num_constrs(m.inner)
        if sense[i] == '='
            # Must be equal to RHS if feasible
            ret[i] = rhs[i]
        else
            # Slack variable is non-negative for <= constraints, and non-positive for >= constraints
            ret[i] = rhs[i] - slack[i]
        end
    end
    return ret
end

function MPB.getreducedcosts(m::GurobiMathProgModel)
    if is_qcp(m.inner) && get_int_param(m.inner, "QCPDual") == 0
        return fill(NaN, num_vars(m.inner))
    else
        return get_dblattrarray(m.inner, "RC", 1, num_vars(m.inner))
    end
end

function MPB.getconstrduals(m::GurobiMathProgModel)
    if is_qcp(m.inner) && get_int_param(m.inner, "QCPDual") == 0
        return fill(NaN, num_constrs(m.inner))
    else
        return get_dblattrarray(m.inner, "Pi", 1, num_constrs(m.inner))
    end
end

function MPB.getquadconstrduals(m::GurobiMathProgModel)
    if is_qcp(m.inner) && get_int_param(m.inner, "QCPDual") == 0
        return fill(NaN, num_qconstrs(m.inner))
    else
        return get_dblattrarray(m.inner, "QCPi", 1, num_qconstrs(m.inner))
    end
end

MPB.getinfeasibilityray(m::GurobiMathProgModel) = -get_dblattrarray(m.inner, "FarkasDual", 1, num_constrs(m.inner)) # note sign is flipped
MPB.getunboundedray(m::GurobiMathProgModel) = get_dblattrarray(m.inner, "UnbdRay", 1, num_vars(m.inner))

MPB.getbasis(m::GurobiMathProgModel) = get_basis(m.inner)

MPB.getrawsolver(m::GurobiMathProgModel) = m.inner

const var_type_map = Dict(
  'C' => :Cont,
  'B' => :Bin,
  'I' => :Int,
  'S' => :SemiCont,
  'N' => :SemiInt
)

const rev_var_type_map = Dict(
  :Cont => 'C',
  :Bin => 'B',
  :Int => 'I',
  :SemiCont => 'S',
  :SemiInt => 'N'
)

function MPB.setvartype!(m::GurobiMathProgModel, vartype::Vector{Symbol})
    # do this to make sure we deal with new columns
    updatemodel!(m)
    nvartype = map(x->rev_var_type_map[x], vartype)
    set_charattrarray!(m.inner, "VType", 1, length(nvartype), nvartype)
    updatemodel!(m) # otherwise getvartype! will return old values
end

function MPB.getvartype(m::GurobiMathProgModel)
    ret = get_charattrarray(m.inner, "VType", 1, num_vars(m.inner))
    map(x->var_type_map[x], ret)
end

function MPB.setwarmstart!(m::GurobiMathProgModel, v)
    for j = 1:length(v)
        if isnan(v[j])
            v[j] = 1e101  # GRB_UNDEFINED
        end
    end
    set_dblattrarray!(m.inner, "Start", 1, num_vars(m.inner), v)
end

MPB.addsos1!(m::GurobiMathProgModel, idx, weight) = add_sos!(m.inner, :SOS1, idx, weight)
MPB.addsos2!(m::GurobiMathProgModel, idx, weight) = add_sos!(m.inner, :SOS2, idx, weight)

# Callbacks


MPB.setlazycallback!(m::GurobiMathProgModel,f) = (m.lazycb = f)
MPB.setcutcallback!(m::GurobiMathProgModel,f) = (m.cutcb = f)
MPB.setheuristiccallback!(m::GurobiMathProgModel,f) = (m.heuristiccb = f)
MPB.setinfocallback!(m::GurobiMathProgModel,f) = (m.infocb = f)

mutable struct GurobiCallbackData <: MPB.MathProgCallbackData
    cbdata::CallbackData
    state::Symbol
    where::Cint
    sol::Vector{Float64}  # Used for heuristic callbacks
#    model::GurobiMathProgModel # not needed?
end

function MPB.cbgetmipsolution(d::GurobiCallbackData)
    @assert d.state == :MIPSol
    return cbget_mipsol_sol(d.cbdata, d.where)
end

function MPB.cbgetmipsolution(d::GurobiCallbackData,output)
    @assert d.state == :MIPSol
    return cbget_mipsol_sol(d.cbdata, d.where, output)
end

function MPB.cbgetlpsolution(d::GurobiCallbackData)
    @assert d.state == :MIPNode
    return cbget_mipnode_rel(d.cbdata, d.where)
end

function MPB.cbgetlpsolution(d::GurobiCallbackData, output)
    @assert d.state == :MIPNode
    return cbget_mipnode_rel(d.cbdata, d.where, output)
end


# TODO: macro for these getters?
function MPB.cbgetobj(d::GurobiCallbackData)
    if d.state == :MIPNode
        return cbget_mipnode_objbst(d.cbdata, d.where)
    elseif d.state == :Intermediate
        return cbget_mip_objbst(d.cbdata, d.where)
    elseif d.state == :MIPSol
        error("Gurobi does not implement cbgetobj when state == MIPSol")
        # https://groups.google.com/forum/#!topic/gurobi/Az_X6Ag-y6k
        # https://github.com/JuliaOpt/MathProgBase.jl/issues/23
        # A complex workaround is possible but hasn't been implemented.
        # return cbget_mipsol_objbst(d.cbdata, d.where)
    else
        error("Unrecognized callback state $(d.state)")
    end
end

function MPB.cbgetbestbound(d::GurobiCallbackData)
    if d.state == :MIPNode
        return cbget_mipnode_objbnd(d.cbdata, d.where)
    elseif d.state == :MIPSol
        return cbget_mipsol_objbnd(d.cbdata, d.where)
    elseif d.state == :Intermediate
        return cbget_mip_objbnd(d.cbdata, d.where)
    else
        error("Unrecognized callback state $(d.state)")
    end
end

function MPB.cbgetexplorednodes(d::GurobiCallbackData)
    if d.state == :MIPNode
        return cbget_mipnode_nodcnt(d.cbdata, d.where)
    elseif d.state == :MIPSol
        return cbget_mipsol_nodcnt(d.cbdata, d.where)
    elseif d.state == :Intermediate
        return cbget_mip_nodcnt(d.cbdata, d.where)
    else
        error("Unrecognized callback state $(d.state)")
    end
end

# returns :MIPNode :MIPSol :Intermediate
MPB.cbgetstate(d::GurobiCallbackData) = d.state

function MPB.cbaddcut!(d::GurobiCallbackData,varidx,varcoef,sense,rhs)
    @assert d.state == :MIPNode
    cbcut(d.cbdata, convert(Vector{Cint}, varidx), float(varcoef), sense, float(rhs))
end

function MPB.cbaddlazy!(d::GurobiCallbackData,varidx,varcoef,sense,rhs)
    @assert d.state == :MIPNode || d.state == :MIPSol
    cblazy(d.cbdata, convert(Vector{Cint}, varidx), float(varcoef), sense, float(rhs))
end

function MPB.cbaddsolution!(d::GurobiCallbackData)
    # Gurobi doesn't support adding solutions on MIPSol.
    # TODO: support this anyway
    @assert d.state == :MIPNode
    cbsolution(d.cbdata, d.sol)
    # "Wipe" solution back to GRB_UNDEFINIED
    for i in 1:length(d.sol)
        d.sol[i] = 1e101  # GRB_UNDEFINED
    end
end

function MPB.cbsetsolutionvalue!(d::GurobiCallbackData,varidx,value)
    d.sol[varidx] = value
end

# breaking abstraction, define our low-level callback to eliminatate
# a level of indirection

function mastercallback(ptr_model::Ptr{Cvoid}, cbdata::Ptr{Cvoid}, where::Cint, userdata::Ptr{Cvoid})

    model = unsafe_pointer_to_objref(userdata)::GurobiMathProgModel
    grbrawcb = CallbackData(cbdata,model.inner)
    if where == CB_MIPSOL
        state = :MIPSol
    elseif where == CB_MIPNODE
        state = :MIPNode
        # skip callback if node is reported to be cut off or infeasible --
        # nothing to do.
        # TODO: users may want this information
        status = cbget_mipnode_status(grbrawcb, where)
        if status != 2
            return convert(Cint,0)
        end
    elseif where == CB_MIP
        state = :Intermediate
    else
        # State with no relevant callbacks
        return convert(Cint,0)
    end

    grbcb = GurobiCallbackData(grbrawcb, state, where, [0.0])
    if model.infocb != nothing
        ret = model.infocb(grbcb)
        if ret == :Exit
            terminate(model.inner)
        end
    end
    if model.cutcb != nothing && state == :MIPNode
        ret = model.cutcb(grbcb)
        if ret == :Exit
            terminate(model.inner)
        end
    end
    if model.heuristiccb != nothing && state == :MIPNode
        grbcb.sol = fill(1e101, MPB.numvar(model))  # GRB_UNDEFINED
        ret = model.heuristiccb(grbcb)
        if ret == :Exit
            terminate(model.inner)
        end
    end
    if model.lazycb != nothing && (state == :MIPSol || state == :MIPNode)
        ret = model.lazycb(grbcb)
        if ret == :Exit
            terminate(model.inner)
        end
    end

    return convert(Cint,0)
end

# User callback function should be of the form:
# callback(cbdata::MPB.MathProgCallbackData)
# return :Exit to indicate an error

function setmathprogcallback!(model::GurobiMathProgModel)
    if Sys.iswindows() && Sys.WORD_SIZE != 64
        error("Callbacks not currently supported on Win32. Use 64-bit Julia with 64-bit Gurobi.")
    end
    grbcallback = @cfunction(mastercallback, Cint, (Ptr{Cvoid}, Ptr{Cvoid}, Cint, Ptr{Cvoid}))
    ret = @grb_ccall(setcallbackfunc, Cint, (Ptr{Cvoid}, Ptr{Cvoid}, Any), model.inner.ptr_model, grbcallback, model)
    if ret != 0
        throw(GurobiError(model.env, ret))
    end
    nothing
end


# QCQP

function MPB.setquadobj!(m::GurobiMathProgModel, rowidx, colidx, quadval)
    delq!(m.inner)
    update_model!(m.inner)
    scaledvals = similar(quadval)
    for i in 1:length(rowidx)
      if rowidx[i] == colidx[i]
        # rescale from matrix format to "terms" format
        scaledvals[i] = quadval[i] / 2
      else
        scaledvals[i] = quadval[i]
      end
    end
    add_qpterms!(m.inner, rowidx, colidx, scaledvals)
    update_model!(m.inner)
end

function MPB.addquadconstr!(m::GurobiMathProgModel, linearidx, linearval, quadrowidx, quadcolidx, quadval, sense, rhs)
    if m.last_op_type == :Var
        updatemodel!(m)
        m.last_op_type = :Con
    end
    add_qconstr!(m.inner, linearidx, linearval, quadrowidx, quadcolidx, quadval, sense, rhs)
end

MPB.getsolvetime(m::GurobiMathProgModel) = get_runtime(m.inner)
MPB.getnodecount(m::GurobiMathProgModel) = get_node_count(m.inner)
MPB.getsimplexiter(m::GurobiMathProgModel) = get_dblattr(m.inner, "IterCount")
