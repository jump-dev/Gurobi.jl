# ==============================================================================
#    Generic Callbacks in Gurobi
# ==============================================================================

"""
    CallbackFunction()

Set a generic Gurobi callback function.

Note: before accessing `MOI.CallbackVariablePrimal`, you must call either
`Gurobi.cbget_mipsol_sol(model, cb_data, cb_where)` or
`Gurobi.cbget_mipsol_rel(model, cb_data, cb_where)`.
"""
struct CallbackFunction <: MOI.AbstractOptimizerAttribute end

function MOI.set(model::Optimizer, ::CallbackFunction, f::Function)
    model.has_generic_callback = true
    set_callback_func!(model.inner, (cb_data, cb_where) -> begin
        model.callback_state = CB_GENERIC
        f(cb_data, cb_where)
        model.callback_state = CB_NONE
    end)
    update_model!(model.inner)
    return
end

"""
    cbget_mipsol_sol(model::Optimizer, cb_data, cb_where)

Load the solution at a CB_MIPSOL node so that it can be accessed using
`MOI.CallbackVariablePrimal`.
"""
function cbget_mipsol_sol(model::Optimizer, cb_data, cb_where)
    resize!(model.callback_variable_primal, length(model.variable_info))
    cbget_mipsol_sol(cb_data, cb_where, model.callback_variable_primal)
    return
end

"""
    cbget_mipsol_rel(model::Optimizer, cb_data, cb_where)

Load the solution at a CB_MIPNODE node so that it can be accessed using
`MOI.CallbackVariablePrimal`.
"""
function cbget_mipsol_rel(model::Optimizer, cb_data, cb_where)
    resize!(model.callback_variable_primal, length(model.variable_info))
    cbget_mipnode_rel(cb_data, cb_where, model.callback_variable_primal)
    return
end

# ==============================================================================
#    MOI callbacks
# ==============================================================================

function default_moi_callback(model::Optimizer)
    return (cb_data, cb_where) -> begin
        if cb_where == CB_MIPSOL
            cbget_mipsol_sol(model, cb_data, cb_where)
            if model.lazy_callback !== nothing
                model.callback_state = CB_LAZY
                model.lazy_callback(cb_data)
            end
        elseif cb_where == CB_MIPNODE
            if cbget_mipnode_status(cb_data, cb_where) != 2
                return  # Solution is something other than optimal.
            end
            cbget_mipsol_rel(model, cb_data, cb_where)
            if model.lazy_callback !== nothing
                model.callback_state = CB_LAZY
                model.lazy_callback(cb_data)
            end
            if model.user_cut_callback !== nothing
                model.callback_state = CB_USER_CUT
                model.user_cut_callback(cb_data)
            end
            if model.heuristic_callback !== nothing
                model.callback_state = CB_HEURISTIC
                model.heuristic_callback(cb_data)
            end
        end
        model.callback_state = CB_NONE
    end
end

function MOI.get(
    model::Optimizer,
    ::MOI.CallbackVariablePrimal{CallbackData},
    x::MOI.VariableIndex
)
    return model.callback_variable_primal[_info(model, x).column]
end

# ==============================================================================
#    MOI.LazyConstraint
# ==============================================================================

function MOI.set(model::Optimizer, ::MOI.LazyConstraintCallback, cb::Function)
    model.lazy_callback = cb
    return
end

function MOI.submit(
    model::Optimizer,
    cb::MOI.LazyConstraint{CallbackData},
    f::MOI.ScalarAffineFunction{Float64},
    s::Union{MOI.LessThan{Float64}, MOI.GreaterThan{Float64}, MOI.EqualTo{Float64}}
)
    if model.callback_state != CB_LAZY && model.callback_state != CB_GENERIC
        error("`MOI.LazyConstraint` can only be called from LazyConstraintCallback.")
    end
    indices, coefficients = _indices_and_coefficients(model, f)
    sense, rhs = _sense_and_rhs(s)
    return cblazy(cb.callback_data, Cint.(indices), coefficients, Char(sense), rhs)
end

# ==============================================================================
#    MOI.UserCutCallback
# ==============================================================================

function MOI.set(model::Optimizer, ::MOI.UserCutCallback, cb::Function)
    model.user_cut_callback = cb
    return
end

function MOI.submit(
    model::Optimizer,
    cb::MOI.UserCut{CallbackData},
    f::MOI.ScalarAffineFunction{Float64},
    s::Union{MOI.LessThan{Float64}, MOI.GreaterThan{Float64}, MOI.EqualTo{Float64}}
)
    if model.callback_state != CB_USER_CUT && model.callback_state != CB_GENERIC
        error("`MOI.UserCut` can only be called from UserCutCallback.")
    end
    indices, coefficients = _indices_and_coefficients(model, f)
    sense, rhs = _sense_and_rhs(s)
    return cbcut(cb.callback_data, Cint.(indices), coefficients, Char(sense), rhs)
end

# ==============================================================================
#    MOI.HeuristicCallback
# ==============================================================================

function MOI.set(model::Optimizer, ::MOI.HeuristicCallback, cb::Function)
    model.heuristic_callback = cb
    return
end

function MOI.submit(
    model::Optimizer,
    cb::MOI.HeuristicSolution{CallbackData},
    variables::Vector{MOI.VariableIndex},
    values::MOI.Vector{Float64}
)
    if model.callback_state != CB_HEURISTIC && model.callback_state != CB_GENERIC
        error("`MOI.HeuristicSolution` can only be called from HeuristicCallback.")
    end
    solution = fill(GRB_UNDEFINED, MOI.get(model, MOI.NumberOfVariables()))
    for (var, value) in zip(variables, values)
        solution[_info(model, var).column] = value
    end
    obj = cbsolution(cb.callback_data, solution)
    return obj < GRB_INFINITY ? MOI.HEURISTIC_SOLUTION_ACCEPTED : MOI.HEURISTIC_SOLUTION_REJECTED
end
