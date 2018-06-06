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
    env::Env
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
LQOI.supported_objectives(s::GurobiOptimizer)  = SUPPORTED_OBJECTIVES

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

function LQOI.get_variable_lowerbound(instance::GurobiOptimizer, column::Int)
    get_dblattrlist(instance.inner, "LB", Cint[column])[1]
end

function LQOI.get_variable_upperbound(instance::GurobiOptimizer, column::Int)
    get_dblattrlist(instance.inner, "UB", Cint[column])[1]
end

function LQOI.get_number_linear_constraints(instance::GurobiOptimizer)
    num_constrs(instance.inner)
end

function LQOI.add_linear_constraints!(instance::GurobiOptimizer,
        A::LQOI.CSRMatrix{Float64}, sense::Vector{Cchar}, rhs::Vector{Float64})
    add_constrs!(instance.inner, A.row_pointers, A.columns, A.coefficients, sense, rhs)
    update_model!(instance.inner)
end

function LQOI.get_rhs(instance::GurobiOptimizer, row::Int)
    get_dblattrlist( instance.inner, "RHS", Cint[row])[1]
end

function LQOI.get_linear_constraint(instance::GurobiOptimizer, row::Int)
    A = get_constrs(instance.inner, row, 1)'
    return A.rowval-1, A.nzval
end

function LQOI.change_matrix_coefficient!(instance::GurobiOptimizer, row::Int, col::Int, coef::Float64)
    chg_coeffs!(instance.inner, row, col, coef)
    update_model!(instance.inner)
end

function LQOI.change_objective_coefficient!(instance::GurobiOptimizer, col::Int, coef::Float64)
    set_dblattrlist!(instance.inner, "Obj", Cint[col], Float64[coef])
    update_model!(instance.inner)
end

function LQOI.change_rhs_coefficient!(instance::GurobiOptimizer, row::Int, coef::Float64)
    set_dblattrlist!(instance.inner, "RHS", Cint[row], Float64[coef])
    update_model!(instance.inner)
end

function LQOI.delete_linear_constraints!(instance::GurobiOptimizer, first_row::Int, last_row::Int)
    del_constrs!(instance.inner, collect(first_row:last_row))
    update_model!(instance.inner)
end

function LQOI.delete_quadratic_constraints!(instance::GurobiOptimizer, first_row::Int, last_row::Int)
    delqconstrs!(instance.inner, collect(first_row:last_row))
    update_model!(instance.inner)
end

function LQOI.change_variable_types!(instance::GurobiOptimizer, columns::Vector{Int}, vtypes::Vector{Cchar})
    set_charattrlist!(instance.inner, "VType", Cint.(columns), vtypes)
    update_model!(instance.inner)
end

function LQOI.change_linear_constraint_sense!(instance::GurobiOptimizer, rows::Vector{Int}, senses::Vector{Cchar})
    set_charattrlist!(instance.inner, "Sense", Cint.(rows), senses)
    update_model!(instance.inner)
end

function LQOI.add_sos_constraint!(instance::GurobiOptimizer, columns::Vector{Int}, weights::Vector{Float64}, sos_type)
    add_sos!(instance.inner, sos_type, columns, weights)
    update_model!(instance.inner)
end

function LQOI.delete_sos!(instance::GurobiOptimizer, first_row::Int, last_row::Int)
    del_sos!(instance.inner, Cint.(first_row:last_row))
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

function LQOI.get_number_quadratic_constraints(instance::GurobiOptimizer)
    num_qconstrs(instance.inner)
end

function scalediagonal!(V, I, J, scale)
    #  LQOI assumes 0.5 x' Q x, but Gurobi requires the list of terms, e.g.,
    #  2x^2 + xy + y^2, so we multiply the diagonal of V by 0.5. We don't
    #  multiply the off-diagonal terms since we assume they are symmetric and we
    #  only need to give one.
    #
    #  We also need to make sure that after adding the constraint we un-scale
    #  the vector because we can't modify user-data.
    for i in 1:length(I)
        if I[i] == J[i]
            V[i] *= scale
        end
    end
end
function LQOI.add_quadratic_constraint!(instance::GurobiOptimizer,
        affine_columns::Vector{Int}, affine_coefficients::Vector{Float64},
        rhs::Float64, sense::Cchar,
        I::Vector{Int}, J::Vector{Int}, V::Vector{Float64})
    @assert length(I) == length(J) == length(V)
    scalediagonal!(V, I, J, 0.5)
    add_qconstr!(instance.inner, affine_columns, affine_coefficients, I, J, V, sense, rhs)
    scalediagonal!(V, I, J, 2.0)
    update_model!(instance.inner)
end

function LQOI.set_quadratic_objective!(instance::GurobiOptimizer, I::Vector{Int}, J::Vector{Int}, V::Vector{Float64})
    @assert length(I) == length(J) == length(V)
    delq!(instance.inner)
    scalediagonal!(V, I, J, 0.5)
    add_qpterms!(instance.inner, I, J, V)
    scalediagonal!(V, I, J, 2.0)
    update_model!(instance.inner)
end

function LQOI.set_linear_objective!(instance::GurobiOptimizer, columns::Vector{Int}, coefficients::Vector{Float64})
    nvars = num_vars(instance.inner)
    obj = zeros(Float64, nvars)
    for (col, coef) in zip(columns, coefficients)
        obj[col] = coef
    end
    set_dblattrarray!(instance.inner, "Obj", 1, num_vars(instance.inner), obj)
    update_model!(instance.inner)
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

function LQOI.get_number_variables(instance::GurobiOptimizer)
    update_model!(instance.inner)
    num_vars(instance.inner)
end

function LQOI.add_variables!(instance::GurobiOptimizer, N::Int)
    add_cvars!(instance.inner, zeros(N))
    update_model!(instance.inner)
end

function LQOI.delete_variables!(instance::GurobiOptimizer, first_col::Int, last_col::Int)
    del_vars!(instance.inner, Cint.(first_col:last_col))
    update_model!(instance.inner)
end

function LQOI.add_mip_starts!(instance::GurobiOptimizer, columns::Vector{Int}, starts::Vector{Float64})
    x = zeros(num_vars(instance.inner))
    for (col, val) in zip(columns, starts)
        x[col] = val
    end
    loadbasis(instance.inner, x)
end

LQOI.solve_mip_problem!(instance::GurobiOptimizer) = LQOI.solve_linear_problem!(instance)

LQOI.solve_quadratic_problem!(instance::GurobiOptimizer) = LQOI.solve_linear_problem!(instance)

function LQOI.solve_linear_problem!(instance::GurobiOptimizer)
    update_model!(instance.inner)
    optimize(instance.inner)
end

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

function LQOI.get_variable_primal_solution!(instance::GurobiOptimizer, result)
    get_dblattrarray!(result, instance.inner, "X", 1)
end

function LQOI.get_linear_primal_solution!(instance::GurobiOptimizer, result)
    get_dblattrarray!(result, instance.inner, "Slack", 1)
    rhs = get_dblattrarray(instance.inner, "RHS", 1, num_constrs(instance.inner))
    result .= rhs - result
end

function LQOI.get_quadratic_primal_solution!(instance::GurobiOptimizer, place)
    get_dblattrarray!(place, instance.inner, "QCSlack", 1)
    rhs = get_dblattrarray(instance.inner, "QCRHS", 1, num_qconstrs(instance.inner))
    place .= rhs - place
end

function LQOI.get_variable_dual_solution!(instance::GurobiOptimizer, place)
    get_dblattrarray!(place, instance.inner, "RC", 1)
end

function LQOI.get_linear_dual_solution!(instance::GurobiOptimizer, place)
    get_dblattrarray!(place, instance.inner, "Pi", 1)
end

function LQOI.get_quadratic_dual_solution!(instance::GurobiOptimizer, place)
    get_dblattrarray!(place, instance.inner, "QCPi", 1)
end

LQOI.get_objective_value(instance::GurobiOptimizer) = get_objval(instance.inner)
LQOI.get_objective_bound(instance::GurobiOptimizer) = get_objval(instance.inner)

function LQOI.get_relative_mip_gap(instance::GurobiOptimizer)
    L = get_objval(instance.inner)
    U = get_objbound(instance.inner)
    return abs(U-L)/U
end

LQOI.get_iteration_count(instance::GurobiOptimizer)    = get_iter_count(instance.inner)
LQOI.get_barrier_iterations(instance::GurobiOptimizer) = get_barrier_iter_count(instance.inner)
LQOI.get_node_count(instance::GurobiOptimizer)         = get_node_count(instance.inner)
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


# ==============================================================================
#    Callbacks in Gurobi
# ==============================================================================
struct CallbackFunction <: MOI.AbstractOptimizerAttribute end
function MOI.set!(m::GurobiOptimizer, ::CallbackFunction, f::Function)
    set_callback_func!(m.inner, f)
    update_model!(m.inner)
end

function cblazy!(cb_data::Gurobi.CallbackData, m::GurobiOptimizer, func::LQOI.Linear, set::S) where S <: Union{LQOI.LE, LQOI.GE, LQOI.EQ}
    columns      = [Cint(LQOI.getcol(m, term.variable_index)) for term in func.terms]
    coefficients = [term.coefficient for term in func.terms]
    sense        = Char(LQOI.backend_type(m, set))
    rhs          = MOI.Utilities.getconstant(set)
    cblazy(cb_data, columns, coefficients, sense, rhs)
end
