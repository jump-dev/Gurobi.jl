# ==============================================================================
#    Generic Callbacks in Gurobi
# ==============================================================================

mutable struct CallbackData
    model::Optimizer
    ptr::Ptr{Cvoid}
    cb_where::Cint
end
Base.cconvert(::Type{Ptr{Cvoid}}, x::CallbackData) = x
Base.unsafe_convert(::Type{Ptr{Cvoid}}, x::CallbackData) = x.ptr::Ptr{Cvoid}
Base.broadcastable(x::CallbackData) = Ref(x)

mutable struct _CallbackUserData
    model::Optimizer
    callback::Function
end
Base.cconvert(::Type{Ptr{Cvoid}}, x::_CallbackUserData) = x
function Base.unsafe_convert(::Type{Ptr{Cvoid}}, x::_CallbackUserData)
    return pointer_from_objref(x)::Ptr{Cvoid}
end

function _gurobi_callback_wrapper(
    p_model::Ptr{Cvoid},
    cb_data::Ptr{Cvoid},
    cb_where::Cint,
    p_user_data::Ptr{Cvoid},
)
    user_data = unsafe_pointer_to_objref(p_user_data)::_CallbackUserData
    try
        user_data.callback(
            CallbackData(user_data.model, cb_data, cb_where),
            cb_where,
        )
    catch ex
        GRBterminate(p_model)
        if !(ex isa InterruptException)
            rethrow(ex)
        end
    end
    return Cint(0)
end

"""
    CallbackFunction()

Set a generic Gurobi callback function.

Callback function should be of the form

    callback(cb_data::CallbackData, cb_where::Cint)

Note: before accessing `MOI.CallbackVariablePrimal`, you must call
`Gurobi.load_callback_variable_primal(cb_data, cb_where)`.
"""
struct CallbackFunction <: MOI.AbstractCallback end

function MOI.set(model::Optimizer, ::CallbackFunction, f::Function)
    grb_callback = @cfunction(
        _gurobi_callback_wrapper,
        Cint,
        (Ptr{Cvoid}, Ptr{Cvoid}, Cint, Ptr{Cvoid})
    )
    user_data = _CallbackUserData(
        model,
        (cb_data, cb_where) -> begin
            model.callback_state = _CB_GENERIC
            f(cb_data, cb_where)
            model.callback_state = _CB_NONE
            return
        end,
    )
    ret = GRBsetcallbackfunc(model, grb_callback, user_data)
    _check_ret(model, ret)
    # We need to keep a reference to the callback function so that it isn't
    # garbage collected.
    model.generic_callback = user_data
    model.has_generic_callback = true
    # Mark the update as necessary and immediately call for the update.
    _require_update(model)
    _update_if_necessary(model)
    return
end
MOI.supports(::Optimizer, ::CallbackFunction) = true

"""
    load_callback_variable_primal(cb_data, cb_where)

Load the solution during a callback so that it can be accessed using
`MOI.CallbackVariablePrimal`.
"""
function load_callback_variable_primal(cb_data::CallbackData, cb_where::Cint)
    resize!(
        cb_data.model.callback_variable_primal,
        length(cb_data.model.variable_info),
    )
    if cb_where == GRB_CB_MIPNODE
        ret = GRBcbget(
            cb_data,
            cb_where,
            GRB_CB_MIPNODE_REL,
            cb_data.model.callback_variable_primal,
        )
        _check_ret(cb_data.model, ret)
    elseif cb_where == GRB_CB_MIPSOL
        ret = GRBcbget(
            cb_data,
            cb_where,
            GRB_CB_MIPSOL_SOL,
            cb_data.model.callback_variable_primal,
        )
        _check_ret(cb_data.model, ret)
    else
        error(
            "`load_callback_variable_primal` can only be called at " *
            "GRB_CB_MIPNODE or GRB_CB_MIPSOL.",
        )
    end
    return
end

# ==============================================================================
#    MOI callbacks
# ==============================================================================

function _default_moi_callback(model::Optimizer)
    return (cb_data, cb_where) -> begin
        if cb_where == GRB_CB_MIPSOL
            load_callback_variable_primal(cb_data, cb_where)
            if model.lazy_callback !== nothing
                model.callback_state = _CB_LAZY
                model.lazy_callback(cb_data)
            end
        elseif cb_where == GRB_CB_MIPNODE
            resultP = Ref{Cint}()
            GRBcbget(cb_data, cb_where, GRB_CB_MIPNODE_STATUS, resultP)
            if resultP[] != GRB_OPTIMAL
                return  # Solution is something other than optimal.
            end
            load_callback_variable_primal(cb_data, cb_where)
            if model.user_cut_callback !== nothing
                model.callback_state = _CB_USER_CUT
                model.user_cut_callback(cb_data)
            end
            if model.heuristic_callback !== nothing
                model.callback_state = _CB_HEURISTIC
                model.heuristic_callback(cb_data)
            end
            if model.lazy_callback !== nothing
                model.callback_state = _CB_LAZY
                model.lazy_callback(cb_data)
            end
        end
        model.callback_state = _CB_NONE
    end
end

function MOI.get(
    model::Optimizer,
    ::MOI.CallbackVariablePrimal{CallbackData},
    x::MOI.VariableIndex,
)
    return model.callback_variable_primal[_info(model, x).column]
end

function MOI.get(::Optimizer, attr::MOI.CallbackNodeStatus{CallbackData})
    if attr.callback_data.cb_where == GRB_CB_MIPSOL
        return MOI.CALLBACK_NODE_STATUS_INTEGER
    elseif attr.callback_data.cb_where == GRB_CB_MIPNODE
        return MOI.CALLBACK_NODE_STATUS_FRACTIONAL
    end
    return MOI.CALLBACK_NODE_STATUS_UNKNOWN
end

# ==============================================================================
#    MOI.LazyConstraint
# ==============================================================================

function MOI.set(model::Optimizer, ::MOI.LazyConstraintCallback, cb::Function)
    MOI.set(model, MOI.RawOptimizerAttribute("LazyConstraints"), 1)
    model.lazy_callback = cb
    return
end
MOI.supports(::Optimizer, ::MOI.LazyConstraintCallback) = true

function MOI.submit(
    model::Optimizer,
    cb::MOI.LazyConstraint{CallbackData},
    f::MOI.ScalarAffineFunction{Float64},
    s::Union{
        MOI.LessThan{Float64},
        MOI.GreaterThan{Float64},
        MOI.EqualTo{Float64},
    },
)
    if model.callback_state == _CB_USER_CUT
        throw(MOI.InvalidCallbackUsage(MOI.UserCutCallback(), cb))
    elseif model.callback_state == _CB_HEURISTIC
        throw(MOI.InvalidCallbackUsage(MOI.HeuristicCallback(), cb))
    elseif !iszero(f.constant)
        throw(
            MOI.ScalarFunctionConstantNotZero{Float64,typeof(f),typeof(s)}(
                f.constant,
            ),
        )
    end
    indices, coefficients = _indices_and_coefficients(model, f)
    sense, rhs = _sense_and_rhs(s)
    ret = GRBcblazy(
        cb.callback_data,
        length(indices),
        indices,
        coefficients,
        sense,
        rhs,
    )
    _check_ret(model, ret)
    return
end
MOI.supports(::Optimizer, ::MOI.LazyConstraint{CallbackData}) = true

# ==============================================================================
#    MOI.UserCutCallback
# ==============================================================================

function MOI.set(model::Optimizer, ::MOI.UserCutCallback, cb::Function)
    model.user_cut_callback = cb
    return
end
MOI.supports(::Optimizer, ::MOI.UserCutCallback) = true

function MOI.submit(
    model::Optimizer,
    cb::MOI.UserCut{CallbackData},
    f::MOI.ScalarAffineFunction{Float64},
    s::Union{
        MOI.LessThan{Float64},
        MOI.GreaterThan{Float64},
        MOI.EqualTo{Float64},
    },
)
    if model.callback_state == _CB_LAZY
        throw(MOI.InvalidCallbackUsage(MOI.LazyConstraintCallback(), cb))
    elseif model.callback_state == _CB_HEURISTIC
        throw(MOI.InvalidCallbackUsage(MOI.HeuristicCallback(), cb))
    elseif !iszero(f.constant)
        throw(
            MOI.ScalarFunctionConstantNotZero{Float64,typeof(f),typeof(s)}(
                f.constant,
            ),
        )
    end
    indices, coefficients = _indices_and_coefficients(model, f)
    sense, rhs = _sense_and_rhs(s)
    ret = GRBcbcut(
        cb.callback_data,
        length(indices),
        indices,
        coefficients,
        sense,
        rhs,
    )
    _check_ret(model, ret)
    return
end
MOI.supports(::Optimizer, ::MOI.UserCut{CallbackData}) = true

# ==============================================================================
#    MOI.HeuristicCallback
# ==============================================================================

function MOI.set(model::Optimizer, ::MOI.HeuristicCallback, cb::Function)
    model.heuristic_callback = cb
    return
end
MOI.supports(::Optimizer, ::MOI.HeuristicCallback) = true

function MOI.submit(
    model::Optimizer,
    cb::MOI.HeuristicSolution{CallbackData},
    variables::Vector{MOI.VariableIndex},
    values::MOI.Vector{Float64},
)
    if model.callback_state == _CB_LAZY
        throw(MOI.InvalidCallbackUsage(MOI.LazyConstraintCallback(), cb))
    elseif model.callback_state == _CB_USER_CUT
        throw(MOI.InvalidCallbackUsage(MOI.UserCutCallback(), cb))
    end
    solution = fill(GRB_UNDEFINED, MOI.get(model, MOI.NumberOfVariables()))
    for (var, value) in zip(variables, values)
        solution[_info(model, var).column] = value
    end
    objP = Ref{Cdouble}()
    ret = GRBcbsolution(cb.callback_data, solution, objP)
    _check_ret(model, ret)
    return objP[] < GRB_INFINITY ? MOI.HEURISTIC_SOLUTION_ACCEPTED :
           MOI.HEURISTIC_SOLUTION_REJECTED
end
MOI.supports(::Optimizer, ::MOI.HeuristicSolution{CallbackData}) = true
