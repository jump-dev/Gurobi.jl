export GurobiOptimizer

using LinQuadOptInterface
const LQOI = LinQuadOptInterface
const MOI  = LQOI.MOI

const SUPPORTED_OBJECTIVES = [
    LQOI.Linear,
    LQOI.Quad
]

const SUPPORTED_CONSTRAINTS = [
    (LQOI.Linear, LQOI.EQ),
    (LQOI.Linear, LQOI.LE),
    (LQOI.Linear, LQOI.GE),
    # (Linear, IV),
    (LQOI.Quad, LQOI.EQ),
    (LQOI.Quad, LQOI.LE),
    (LQOI.Quad, LQOI.GE),
    (LQOI.SinVar, LQOI.EQ),
    (LQOI.SinVar, LQOI.LE),
    (LQOI.SinVar, LQOI.GE),
    (LQOI.SinVar, LQOI.IV),
    (LQOI.SinVar, MOI.ZeroOne),
    (LQOI.SinVar, MOI.Integer),
    (LQOI.VecVar, LQOI.SOS1),
    (LQOI.VecVar, LQOI.SOS2),
    (LQOI.SinVar, MOI.Semicontinuous{Float64}),
    (LQOI.SinVar, MOI.Semiinteger{Float64}),
    (LQOI.VecVar, MOI.Nonnegatives),
    (LQOI.VecVar, MOI.Nonpositives),
    (LQOI.VecVar, MOI.Zeros),
    (LQOI.VecLin, MOI.Nonnegatives),
    (LQOI.VecLin, MOI.Nonpositives),
    (LQOI.VecLin, MOI.Zeros)
]

mutable struct GurobiOptimizer <: LQOI.LinQuadOptimizer
    LQOI.@LinQuadOptimizerBase
    env
    params::Dict{String,Any}
    GurobiOptimizer(::Void) = new()
end

LQOI.LinearQuadraticModel(::Type{GurobiOptimizer},env) = Model(env::Env,"defaultname")

function GurobiOptimizer(;kwargs...)

    env = Env()
    m = GurobiOptimizer(nothing)
    m.env = env
    m.params = Dict{String,Any}()
    MOI.empty!(m)
    for (name,value) in kwargs
        m.params[string(name)] = value
        setparam!(m.inner, string(name), value)
    end
    return m
end

function MOI.empty!(m::GurobiOptimizer)
    MOI.empty!(m,m.env)
    for (name,value) in m.params
        setparam!(m.inner, name, value)
    end
end

LQOI.supported_constraints(s::GurobiOptimizer) = SUPPORTED_CONSTRAINTS
LQOI.supported_objectives(s::GurobiOptimizer) = SUPPORTED_OBJECTIVES

cintvec(v::Vector) = convert(Vector{Int32}, v)

LQOI.backend_type(m::GurobiOptimizer, ::MOI.EqualTo{Float64})     = Cchar('=')
LQOI.backend_type(m::GurobiOptimizer, ::MOI.LessThan{Float64})    = Cchar('<')
LQOI.backend_type(m::GurobiOptimizer, ::MOI.GreaterThan{Float64}) = Cchar('>')
LQOI.backend_type(m::GurobiOptimizer, ::MOI.Zeros)                = Cchar('=')
LQOI.backend_type(m::GurobiOptimizer, ::MOI.Nonpositives)         = Cchar('<')
LQOI.backend_type(m::GurobiOptimizer, ::MOI.Nonnegatives)         = Cchar('>')

# TODO - improve single type
function LQOI.change_variable_bounds!(instance::GurobiOptimizer, colvec, valvec, sensevec)
    lb_len = count(x->x==Cchar('L'), sensevec)
    LB_val = Array{Float64}(0)
    sizehint!(LB_val, lb_len)
    LB_col = Array{Cint}(0)
    sizehint!(LB_col, lb_len)

    ub_len = count(x->x==Cchar('U'), sensevec)
    UB_val = Array{Float64}(0)
    sizehint!(UB_val, ub_len)
    UB_col = Array{Cint}(0)
    sizehint!(UB_col, ub_len)

    for i in eachindex(valvec)
        if sensevec[i] == Cchar('L')
            push!(LB_col, colvec[i])
            push!(LB_val, valvec[i])
        elseif sensevec[i] == Cchar('U')
            push!(UB_col, colvec[i])
            push!(UB_val, valvec[i])
        end
    end

    if lb_len > 0
        set_dblattrlist!(instance.inner, "LB", LB_col, LB_val)
    end

    if ub_len > 0
        set_dblattrlist!(instance.inner, "UB", UB_col, UB_val)
    end

    update_model!(instance.inner)
    nothing
end

function LQOI.get_variable_lowerbound(instance::GurobiOptimizer, col)
    # TODO(odow): is this needed? update_model!(instance.inner)
    get_dblattrlist(instance.inner, "LB", ivec(col))[1]
end

function LQOI.get_variable_upperbound(instance::GurobiOptimizer, col)
    get_dblattrlist(instance.inner, "UB", ivec(col))[1]
end

function LQOI.get_number_linear_constraints(instance::GurobiOptimizer)
    num_constrs(instance.inner)
end

function LQOI.add_linear_constraints!(instance::GurobiOptimizer, rowvec, colvec, coefvec, sensevec, rhsvec)
    add_constrs!(instance.inner, rowvec, colvec, coefvec, sensevec, rhsvec)
    update_model!(instance.inner)
end

function LQOI.get_rhs(instance::GurobiOptimizer, row)
    get_dblattrlist( instance.inner, "RHS", ivec(row))[1]
end

function LQOI.get_linear_constraint(instance::GurobiOptimizer, idx)
    A = get_constrs(instance.inner, idx, 1)'
    return A.rowval-1, A.nzval
end

# TODO SPLIT THIS ONE
function LQOI.change_coefficient!(instance::GurobiOptimizer, row, col, coef)
    if row == 0
        set_dblattrlist!(instance.inner, "Obj", Cint[col], Float64[coef])
    elseif col == 0
        set_dblattrlist!(instance.inner, "RHS", Cint[row], Float64[coef])
    else
        chg_coeffs!(instance.inner, row, col, coef)
        #TODO fix this function in gurobi
    end
end

function LQOI.delete_linear_constraints!(instance::GurobiOptimizer, rowbeg, rowend)
    del_constrs!(instance.inner, cintvec(collect(rowbeg:rowend)))
    update_model!(instance.inner)
end

function LQOI.delete_quadratic_constraints!(instance::GurobiOptimizer, i, j)
    delqconstrs!(instance.inner, collect(i:j))
    update_model!(instance.inner)
end

# TODO fix types
LQOI.change_variable_types!(instance::GurobiOptimizer, colvec, typevec) = set_charattrlist!(instance.inner, "VType", ivec(colvec), cvec(typevec))

# TODO fix types
LQOI.change_linear_constraint_sense!(instance::GurobiOptimizer, rowvec, sensevec) = set_charattrlist!(instance.inner, "Sense", ivec(rowvec), cvec(sensevec))

function LQOI.add_sos_constraint!(instance::GurobiOptimizer, colvec, valvec, typ)
    add_sos!(instance.inner, typ, colvec, valvec)
    update_model!(instance.inner)
end

function LQOI.delete_sos!(instance::GurobiOptimizer, idx1, idx2)
    del_sos!(instance.inner, cintvec(collect(idx1:idx2)))
    update_model!(instance.inner)
end

# TODO improve getting processes
function LQOI.get_sos_constraint(instance::GurobiOptimizer, idx)
    A, types = get_sos_matrix(instance.inner)
    line = A[idx,:] #sparse vec
    cols = line.nzind
    vals = line.nzval
    typ = types[idx] == Cint(1) ? :SOS1 : :SOS2
    return cols, vals, typ
end

LQOI.get_number_quadratic_constraints(instance::GurobiOptimizer) = num_qconstrs(instance.inner)

#   NOTE:
# LQOI assumes 0.5 x' Q x, but Gurobi requires x' Q x so we multiply V by 0.5
function LQOI.add_quadratic_constraint!(instance::GurobiOptimizer, cols,coefs,rhs,sense, I,J,V)
    add_qconstr!(instance.inner, cols, coefs, I, J, 0.5 * V, sense, rhs)
    update_model!(instance.inner)
end

# LQOI.change_range_value!(instance::GurobiOptimizer, rows, vals) = chg_rhsrange!(instance.inner, cintvec(rows), -vals)

function LQOI.set_quadratic_objective!(instance::GurobiOptimizer, I, J, V)
    delq!(instance.inner)
    for i in eachindex(V)
        if I[i] == J[i]
            V[i] /= 2
        end
    end
    add_qpterms!(instance.inner, I, J, V)
    update_model!(instance.inner)
    return nothing
end

function LQOI.set_linear_objective!(instance::GurobiOptimizer, colvec, coefvec)
    nvars = num_vars(instance.inner)
    obj = zeros(Float64, nvars)

    for i in eachindex(colvec)
        obj[colvec[i]] = coefvec[i]
    end

    set_dblattrarray!(instance.inner, "Obj", 1, num_vars(instance.inner), obj)
    update_model!(instance.inner)
    nothing
end

function LQOI.change_objective_sense!(instance::GurobiOptimizer, symbol)
    if symbol == :min
        set_sense!(instance.inner, :minimize)
    else
        set_sense!(instance.inner, :maximize)
    end
    update_model!(instance.inner)
end

function LQOI.get_linear_objective!(instance::GurobiOptimizer, x)
    copy!(x, get_dblattrarray(instance.inner, "Obj", 1, num_vars(instance.inner)))
end

function LQOI.get_objectivesense(instance::GurobiOptimizer)
    s = model_sense(instance.inner)
    if s == :maximize
        return MOI.MaxSense
    else
        return MOI.MinSense
    end
end

LQOI.get_number_variables(instance::GurobiOptimizer) = (update_model!(instance.inner); num_vars(instance.inner))

LQOI.add_variables!(instance::GurobiOptimizer, int) = (add_cvars!(instance.inner, zeros(int));update_model!(instance.inner))

# TODO(odow): is this implemented correctly?
LQOI.delete_variables!(instance::GurobiOptimizer, col, col2) = (del_vars!(instance.inner, col);update_model!(instance.inner))

function LQOI.add_mip_starts!(instance::GurobiOptimizer, colvec::Vector{Int}, valvec::Vector)
    x = zeros(num_vars(instance.inner))
    for (col, val) in zip(colvec, valvec)
        x[col] = val
    end
    loadbasis(instance.inner, x)
end

LQOI.solve_mip_problem!(instance::GurobiOptimizer) = LQOI.solve_linear_problem!(instance)

LQOI.solve_quadratic_problem!(instance::GurobiOptimizer) = LQOI.solve_linear_problem!(instance)

LQOI.solve_linear_problem!(instance::GurobiOptimizer) = (update_model!(instance.inner);optimize(instance.inner))

function LQOI.get_termination_status(instance::GurobiOptimizer)
    stat = get_status(instance.inner)
    if stat == :loaded
        return MOI.OtherError
    elseif stat == :optimal
        return MOI.Success
    elseif stat == :infeasible
        if hasdualray(instance)
            return MOI.Success
        else
            return MOI.InfeasibleNoResult
        end
    elseif stat == :inf_or_unbd
        return MOI.InfeasibleOrUnbounded
    elseif stat == :unbounded
        if hasprimalray(instance)
            return MOI.Success
        else
            return MOI.UnboundedNoResult
        end
    elseif stat == :cutoff
        return MOI.ObjectiveLimit
    elseif stat == :iteration_limit
        return MOI.IterationLimit
    elseif stat == :node_limit
        return MOI.NodeLimit
    elseif stat == :time_limit
        return MOI.TimeLimit
    elseif stat == :solution_limit
        return MOI.SolutionLimit
    elseif stat == :interrupted
        return MOI.Interrupted
    elseif stat == :numeric
        return MOI.NumericalError
    elseif stat == :suboptimal
        return MOI.OtherLimit
    elseif stat == :inprogress
        return MOI.OtherError
    elseif stat == :user_obj_limit
        return MOI.ObjectiveLimit
    end
    return MOI.OtherError
end

function LQOI.get_primal_status(instance::GurobiOptimizer)

    stat = get_status(instance.inner)

    if stat == :optimal
        return MOI.FeasiblePoint
    elseif stat == :solution_limit
        return MOI.FeasiblePoint
    elseif stat in [:inf_or_unbd, :unbounded] && hasprimalray(instance)
        return MOI.InfeasibilityCertificate
    elseif stat == :suboptimal
        return MOI.FeasiblePoint
    else
        return MOI.UnknownResultStatus
    end
end

function LQOI.get_dual_status(instance::GurobiOptimizer)
    stat = get_status(instance.inner)

    if is_mip(instance.inner) || is_qcp(instance.inner)
        return MOI.UnknownResultStatus
    else
        if stat == :optimal
            return MOI.FeasiblePoint
        elseif stat == :solution_limit
            return MOI.FeasiblePoint
        elseif stat in [:inf_or_unbd, :infeasible] && hasdualray(instance)
            return MOI.InfeasibilityCertificate
        elseif stat == :suboptimal
            return MOI.FeasiblePoint
        else
            return MOI.UnknownResultStatus
        end
    end
end

LQOI.get_variable_primal_solution!(instance::GurobiOptimizer, place) = get_dblattrarray!(place, instance.inner, "X", 1)

function LQOI.get_linear_primal_solution!(instance::GurobiOptimizer, place)
    get_dblattrarray!(place, instance.inner, "Slack", 1)
    rhs = get_dblattrarray(instance.inner, "RHS", 1, num_constrs(instance.inner))
    for i in eachindex(place)
        place[i] = -place[i]+rhs[i]
    end
    nothing
end

function LQOI.get_quadratic_primal_solution!(instance::GurobiOptimizer, place)
    get_dblattrarray!(place, instance.inner, "QCSlack", 1)
    rhs = get_dblattrarray(instance.inner, "QCRHS", 1, num_qconstrs(instance.inner))
    for i in eachindex(place)
        place[i] = -place[i]+rhs[i]
    end
    nothing
end

LQOI.get_variable_dual_solution!(instance::GurobiOptimizer, place) = get_dblattrarray!(place, instance.inner, "RC", 1)

LQOI.get_linear_dual_solution!(instance::GurobiOptimizer, place) = get_dblattrarray!(place, instance.inner, "Pi", 1)

LQOI.get_quadratic_dual_solution!(instance::GurobiOptimizer, place) = get_dblattrarray!(place, instance.inner, "QCPi", 1)

LQOI.get_objective_value(instance::GurobiOptimizer) = get_objval(instance.inner)

LQOI.get_objective_bound(instance::GurobiOptimizer) = get_objval(instance.inner)

function LQOI.get_relative_mip_gap(instance::GurobiOptimizer)
    L = get_objval(instance.inner)
    U = get_objbound(instance.inner)
    return abs(U-L)/U
end

LQOI.get_iteration_count(instance::GurobiOptimizer)  = get_iter_count(instance.inner)

LQOI.get_barrier_iterations(instance::GurobiOptimizer) = get_barrier_iter_count(instance.inner)

LQOI.get_node_count(instance::GurobiOptimizer) = get_node_count(instance.inner)

LQOI.get_farkas_dual!(instance::GurobiOptimizer, place) = get_dblattrarray!(place, instance.inner, "FarkasDual", 1)

function hasdualray(instance::GurobiOptimizer)
    try
        get_dblattrarray(instance.inner, "FarkasDual", 1, num_constrs(instance.inner))
        return true
    catch
        return false
    end
end

LQOI.get_unbounded_ray!(instance::GurobiOptimizer, place) = get_dblattrarray!(place, instance.inner, "UnbdRay", 1)

function hasprimalray(instance::GurobiOptimizer)
    try
        get_dblattrarray(instance.inner, "UnbdRay", 1, num_vars(instance.inner))
        return true
    catch
        return false
    end
end

MOI.free!(m::GurobiOptimizer) = free_model(m.inner)
