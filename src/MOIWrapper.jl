using LinQuadOptInterface

const LQOI = LinQuadOptInterface
const MOI  = LQOI.MOI

const SUPPORTED_OBJECTIVES = [
    LQOI.SinVar,
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

mutable struct Optimizer <: LQOI.LinQuadOptimizer
    LQOI.@LinQuadOptimizerBase(Model)
    env::Env
    params::Dict{String,Any}
    Optimizer(::Nothing) = new()
end

LQOI.LinearQuadraticModel(::Type{Optimizer}, env) = Model(env::Env,"defaultname")

"""
    Optimizer(;kwargs...)

Create a new Optimizer object.

Note that we set the parameter `InfUnbdInfo` to `1` rather than the default of
`0` so that we can query infeasibility certificates. Users are, however, free to
overide this as follows `Gurobi.Optimizer(InfUndbInfo=0)`.
"""
function Optimizer(;kwargs...)
    model = Optimizer(nothing)
    model.env = Env()
    model.params = Dict{String,Any}()
    MOI.empty!(model)
    for (name, value) in kwargs
        model.params[string(name)] = value
        setparam!(model.inner, string(name), value)
    end
    return model
end

function MOI.empty!(model::Optimizer)
    MOI.empty!(model, model.env)
    setparam!(model.inner, "InfUnbdInfo", 1)
    for (name, value) in model.params
        setparam!(model.inner, name, value)
    end
end

LQOI.supported_constraints(::Optimizer) = SUPPORTED_CONSTRAINTS
LQOI.supported_objectives(::Optimizer)  = SUPPORTED_OBJECTIVES

LQOI.backend_type(::Optimizer, ::MOI.EqualTo{Float64})     = Cchar('=')
LQOI.backend_type(::Optimizer, ::MOI.LessThan{Float64})    = Cchar('<')
LQOI.backend_type(::Optimizer, ::MOI.GreaterThan{Float64}) = Cchar('>')
LQOI.backend_type(::Optimizer, ::MOI.Zeros)                = Cchar('=')
LQOI.backend_type(::Optimizer, ::MOI.Nonpositives)         = Cchar('<')
LQOI.backend_type(::Optimizer, ::MOI.Nonnegatives)         = Cchar('>')

function LQOI.change_variable_bounds!(model::Optimizer,
          columns::Vector{Int}, new_bounds::Vector{Float64},
          senses::Vector{Cchar})
    number_lower_bounds = count(x->x==Cchar('L'), senses)
    lower_cols   = fill(0, number_lower_bounds)
    lower_values = fill(0.0, number_lower_bounds)
    number_upper_bounds = count(x->x==Cchar('U'), senses)
    upper_cols   = fill(0, number_upper_bounds)
    upper_values = fill(0.0, number_upper_bounds)
    lower_index = 1
    upper_index = 1
    for (column, bound, sense) in zip(columns, new_bounds, senses)
        if sense == Cchar('L')
            lower_cols[lower_index]   = column
            lower_values[lower_index] = bound
            lower_index += 1
        elseif sense == Cchar('U')
            upper_cols[upper_index]   = column
            upper_values[upper_index] = bound
            upper_index += 1
        end
    end
    if number_lower_bounds > 0
        set_dblattrlist!(model.inner, "LB", lower_cols, lower_values)
    end
    if number_upper_bounds > 0
        set_dblattrlist!(model.inner, "UB", upper_cols, upper_values)
    end
    update_model!(model.inner)
    return
end

function LQOI.get_variable_lowerbound(model::Optimizer, column::Int)
    return get_dblattrelement(model.inner, "LB", column)
end

function LQOI.get_variable_upperbound(model::Optimizer, column::Int)
    return get_dblattrelement(model.inner, "UB", column)
end

function LQOI.get_number_linear_constraints(model::Optimizer)
    return num_constrs(model.inner)
end

function LQOI.add_linear_constraints!(model::Optimizer,
        A::LQOI.CSRMatrix{Float64}, sense::Vector{Cchar}, rhs::Vector{Float64})
    add_constrs!(model.inner, A.row_pointers, A.columns, A.coefficients, sense, rhs)
    update_model!(model.inner)
    return
end

function LQOI.get_rhs(model::Optimizer, row::Int)
    return get_dblattrelement(model.inner, "RHS", row)
end

function LQOI.get_linear_constraint(model::Optimizer, row::Int)
    A = sparse(get_constrs(model.inner, row, 1)')
    # note: we return 1-index columns
    return A.rowval, A.nzval
end

function LQOI.change_matrix_coefficient!(model::Optimizer, row::Int, col::Int, coef::Float64)
    chg_coeffs!(model.inner, row, col, coef)
    update_model!(model.inner)
    return
end

function LQOI.change_objective_coefficient!(model::Optimizer, col::Int, coef::Float64)
    set_dblattrelement!(model.inner, "Obj", col, coef)
    update_model!(model.inner)
    return
end

function LQOI.change_rhs_coefficient!(model::Optimizer, row::Int, coef::Float64)
    set_dblattrelement!(model.inner, "RHS", row, coef)
    update_model!(model.inner)
    return
end

function LQOI.delete_linear_constraints!(model::Optimizer, first_row::Int, last_row::Int)
    del_constrs!(model.inner, collect(first_row:last_row))
    update_model!(model.inner)
    return
end

function LQOI.delete_quadratic_constraints!(model::Optimizer, first_row::Int, last_row::Int)
    delqconstrs!(model.inner, collect(first_row:last_row))
    update_model!(model.inner)
    return
end

function LQOI.change_variable_types!(model::Optimizer, columns::Vector{Int}, vtypes::Vector{Cchar})
    set_charattrlist!(model.inner, "VType", Cint.(columns), vtypes)
    update_model!(model.inner)
    return
end

function LQOI.change_linear_constraint_sense!(model::Optimizer, rows::Vector{Int}, senses::Vector{Cchar})
    set_charattrlist!(model.inner, "Sense", Cint.(rows), senses)
    update_model!(model.inner)
    return
end

function LQOI.add_sos_constraint!(model::Optimizer, columns::Vector{Int}, weights::Vector{Float64}, sos_type)
    add_sos!(model.inner, sos_type, columns, weights)
    update_model!(model.inner)
    return
end

function LQOI.delete_sos!(model::Optimizer, first_row::Int, last_row::Int)
    del_sos!(model.inner, Cint.(first_row:last_row))
    update_model!(model.inner)
    return
end

# TODO improve getting processes
function LQOI.get_sos_constraint(model::Optimizer, idx)
    A, types = get_sos_matrix(model.inner)
    line = A[idx,:] #sparse vec
    cols = line.nzind
    vals = line.nzval
    typ = types[idx] == Cint(1) ? :SOS1 : :SOS2
    return cols, vals, typ
end

function LQOI.get_number_quadratic_constraints(model::Optimizer)
    return num_qconstrs(model.inner)
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
    return
end
function LQOI.add_quadratic_constraint!(model::Optimizer,
        affine_columns::Vector{Int}, affine_coefficients::Vector{Float64},
        rhs::Float64, sense::Cchar,
        I::Vector{Int}, J::Vector{Int}, V::Vector{Float64})
    @assert length(I) == length(J) == length(V)
    scalediagonal!(V, I, J, 0.5)
    add_qconstr!(model.inner, affine_columns, affine_coefficients, I, J, V, sense, rhs)
    scalediagonal!(V, I, J, 2.0)
    update_model!(model.inner)
    return
end

function LQOI.get_quadratic_constraint(model::Optimizer, row::Int)
    affine_cols, affine_coefficients, I, J, V = getqconstr(model.inner, row)
    # note: we return 1-index columns here
    affine_cols .+= 1
    I .+= 1
    J .+= 1
    return affine_cols, affine_coefficients, sparse(I, J, V)
end

function LQOI.get_quadratic_rhs(model::Optimizer, row::Int)
    return get_dblattrelement(model.inner, "QCRHS", row)
end

function LQOI.set_quadratic_objective!(model::Optimizer, I::Vector{Int}, J::Vector{Int}, V::Vector{Float64})
    @assert length(I) == length(J) == length(V)
    delq!(model.inner)
    scalediagonal!(V, I, J, 0.5)
    add_qpterms!(model.inner, I, J, V)
    scalediagonal!(V, I, J, 2.0)
    update_model!(model.inner)
    return
end

function LQOI.set_linear_objective!(model::Optimizer, columns::Vector{Int}, coefficients::Vector{Float64})
    nvars = num_vars(model.inner)
    obj = zeros(Float64, nvars)
    for (col, coef) in zip(columns, coefficients)
        obj[col] += coef
    end
    set_dblattrarray!(model.inner, "Obj", 1, num_vars(model.inner), obj)
    update_model!(model.inner)
    return
end

function LQOI.set_constant_objective!(model::Optimizer, value::Real)
    set_dblattr!(model.inner, "ObjCon", value)
    if num_vars(model.inner) > 0
        # Work-around for https://github.com/JuliaOpt/LinQuadOptInterface.jl/pull/44#issuecomment-409373755
        set_dblattrarray!(model.inner, "Obj", 1, 1,
            get_dblattrarray(model.inner, "Obj", 1, 1))
    end
    update_model!(model.inner)
    return
end

function LQOI.change_objective_sense!(model::Optimizer, sense::Symbol)
    if sense == :min
        set_sense!(model.inner, :minimize)
    elseif sense == :max
        set_sense!(model.inner, :maximize)
    else
        error("Invalid objective sense: $(sense)")
    end
    update_model!(model.inner)
    return
end

function LQOI.get_linear_objective!(model::Optimizer, dest)
    get_dblattrarray!(dest, model.inner, "Obj", 1)
    return dest
end

function LQOI.get_constant_objective(model::Optimizer)
    return get_dblattr(model.inner, "ObjCon")
end

function LQOI.get_objectivesense(model::Optimizer)
    sense = model_sense(model.inner)
    if sense == :maximize
        return MOI.MaxSense
    elseif sense == :minimize
        return MOI.MinSense
    else
        error("Invalid objective sense: $(sense)")
    end
end

function LQOI.get_number_variables(model::Optimizer)
    return num_vars(model.inner)
end

function LQOI.add_variables!(model::Optimizer, N::Int)
    add_cvars!(model.inner, zeros(N))
    update_model!(model.inner)
    return
end

function LQOI.delete_variables!(model::Optimizer, first_col::Int, last_col::Int)
    del_vars!(model.inner, Cint.(first_col:last_col))
    update_model!(model.inner)
    return
end

function LQOI.add_mip_starts!(model::Optimizer, columns::Vector{Int}, starts::Vector{Float64})
    x = zeros(num_vars(model.inner))
    for (col, val) in zip(columns, starts)
        x[col] = val
    end
    loadbasis(model.inner, x)
    update_model!(model.inner)
    return
end

LQOI.solve_mip_problem!(model::Optimizer) = LQOI.solve_linear_problem!(model)

LQOI.solve_quadratic_problem!(model::Optimizer) = LQOI.solve_linear_problem!(model)

function LQOI.solve_linear_problem!(model::Optimizer)
    update_model!(model.inner)
    optimize(model.inner)
    return
end

function LQOI.get_termination_status(model::Optimizer)
    stat = get_status(model.inner)
    if stat == :loaded
        return MOI.OtherError
    elseif stat == :optimal
        return MOI.Success
    elseif stat == :infeasible
        if has_dual_ray(model)
            return MOI.Success
        else
            return MOI.InfeasibleNoResult
        end
    elseif stat == :inf_or_unbd
        return MOI.InfeasibleOrUnbounded
    elseif stat == :unbounded
        if has_primal_ray(model)
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

function LQOI.get_primal_status(model::Optimizer)
    stat = get_status(model.inner)
    if stat == :optimal
        return MOI.FeasiblePoint
    elseif stat == :solution_limit
        return MOI.FeasiblePoint
    elseif (stat == :inf_or_unbd || stat == :unbounded) && has_primal_ray(model)
        return MOI.InfeasibilityCertificate
    elseif stat == :suboptimal
        return MOI.FeasiblePoint
    elseif is_mip(model.inner) && get_sol_count(model.inner) > 0
        return MOI.FeasiblePoint
    else
        return MOI.NoSolution
    end
end

function LQOI.get_dual_status(model::Optimizer)
    stat = get_status(model.inner)
    if is_mip(model.inner) || is_qcp(model.inner)
        return MOI.UnknownResultStatus
    else
        if stat == :optimal
            return MOI.FeasiblePoint
        elseif stat == :solution_limit
            return MOI.FeasiblePoint
        elseif (stat == :inf_or_unbd || stat == :infeasible) && has_dual_ray(model)
            return MOI.InfeasibilityCertificate
        elseif stat == :suboptimal
            return MOI.FeasiblePoint
        end
    end
    return MOI.UnknownResultStatus
end

function LQOI.get_variable_primal_solution!(model::Optimizer, dest)
    get_dblattrarray!(dest, model.inner, "X", 1)
    return
end

function LQOI.get_linear_primal_solution!(model::Optimizer, dest)
    get_dblattrarray!(dest, model.inner, "RHS", 1)
    dest .-= get_dblattrarray(model.inner, "Slack", 1, num_constrs(model.inner))
    return
end

function LQOI.get_quadratic_primal_solution!(model::Optimizer, dest)
    get_dblattrarray!(dest, model.inner, "QCRHS", 1)
    dest .-= get_dblattrarray(model.inner, "QCSlack", 1, num_qconstrs(model.inner))
    return
end

function LQOI.get_variable_dual_solution!(model::Optimizer, dest)
    get_dblattrarray!(dest, model.inner, "RC", 1)
    return
end

function LQOI.get_linear_dual_solution!(model::Optimizer, dest)
    get_dblattrarray!(dest, model.inner, "Pi", 1)
    return
end

function LQOI.get_quadratic_dual_solution!(model::Optimizer, dest)
    get_dblattrarray!(dest, model.inner, "QCPi", 1)
    return
end

LQOI.get_objective_value(model::Optimizer) = get_objval(model.inner)

function LQOI.get_objective_bound(model::Optimizer)
    return get_objbound(model.inner)
end

function LQOI.get_relative_mip_gap(model::Optimizer)
    value = LQOI.get_objective_value(model)
    bound = LQOI.get_objective_bound(model)
    return abs(value - bound) / abs(bound)
end

function LQOI.get_iteration_count(instance::Optimizer)
    return get_iter_count(instance.inner)
end

function LQOI.get_barrier_iterations(instance::Optimizer)
    return get_barrier_iter_count(instance.inner)
end

function LQOI.get_node_count(instance::Optimizer)
    return get_node_count(instance.inner)
end

function LQOI.get_farkas_dual!(instance::Optimizer, dest)
    get_dblattrarray!(dest, instance.inner, "FarkasDual", 1)
    dest .*= -1.0
    return
end

function has_dual_ray(model::Optimizer)
    try
        # Note: for performance reasons, we try to get 1 element because for
        # some versions of Gurobi, we cannot query 0 elements without error.
        Gurobi.get_dblattrarray(model.inner, "FarkasDual", 1, 1)
        return true
    catch ex
        if isa(ex, Gurobi.GurobiError)
            return false
        else
            rethrow(ex)
        end
    end
end

function LQOI.get_unbounded_ray!(model::Optimizer, dest)
    get_dblattrarray!(dest, model.inner, "UnbdRay", 1)
    return
end

function has_primal_ray(model::Optimizer)
    try
        # Note: for performance reasons, we try to get 1 element because for
        # some versions of Gurobi, we cannot query 0 elements without error.
        Gurobi.get_dblattrarray(model.inner, "UnbdRay", 1, 1)
        return true
    catch ex
        if isa(ex, Gurobi.GurobiError)
            return false
        else
            rethrow(ex)
        end
    end
end

# ==============================================================================
#    Callbacks in Gurobi
# ==============================================================================
struct CallbackFunction <: MOI.AbstractOptimizerAttribute end
function MOI.set(model::Optimizer, ::CallbackFunction, f::Function)
    set_callback_func!(model.inner, f)
    update_model!(model.inner)
    return
end

"""
    loadcbsolution!(m::Optimizer, cb_data::GurobiCallbackData, cb_where::Int)

Load the variable primal solution in a callback.

This can only be called in a callback from `CB_MIPSOL`. After it is called, you
can access the `VariablePrimal` attribute as usual.
"""
function loadcbsolution!(model::Optimizer, cb_data::CallbackData, cb_where::Cint)
    if cb_where != CB_MIPSOL
        error("loadcbsolution! can only be called from CB_MIPSOL.")
    end
    Gurobi.cbget_mipsol_sol(cb_data, cb_where, model.variable_primal_solution)
    return
end

"""
    cblazy!(cb_data::Gurobi.CallbackData, m::Optimizer, func::LQOI.Linear, set::S) where S <: Union{LQOI.LE, LQOI.GE, LQOI.EQ}

Add a lazy cut to the model `m`.

You must have the option `LazyConstraints` set  via `Optimizer(LazyConstraint=1)`.
This can only be called in a callback from `CB_MIPSOL`.
"""
function cblazy!(cb_data::CallbackData, model::Optimizer,
        func::LQOI.Linear, set::S) where S <: Union{LQOI.LE, LQOI.GE, LQOI.EQ}
    columns = [
        Cint(LQOI.get_column(model, term.variable_index)) for term in func.terms]
    coefficients = [term.coefficient for term in func.terms]
    sense = Char(LQOI.backend_type(model, set))
    rhs = MOI.Utilities.getconstant(set)
    return cblazy(cb_data, columns, coefficients, sense, rhs)
end
