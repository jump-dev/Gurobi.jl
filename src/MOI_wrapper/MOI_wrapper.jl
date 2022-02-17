import MathOptInterface

const MOI = MathOptInterface
const CleverDicts = MOI.Utilities.CleverDicts

@enum(
    _BoundType,
    _NONE,
    _LESS_THAN,
    _GREATER_THAN,
    _LESS_AND_GREATER_THAN,
    _INTERVAL,
    _EQUAL_TO
)
@enum(_ObjectiveType, _SINGLE_VARIABLE, _SCALAR_AFFINE, _SCALAR_QUADRATIC)
@enum(
    _CallbackState,
    _CB_NONE,
    _CB_GENERIC,
    _CB_LAZY,
    _CB_USER_CUT,
    _CB_HEURISTIC
)

const _SCALAR_SETS = Union{
    MOI.GreaterThan{Float64},
    MOI.LessThan{Float64},
    MOI.EqualTo{Float64},
    MOI.Interval{Float64},
}

# Union used by many methods because interval constraints are not supported.
const _SUPPORTED_SCALAR_SETS =
    Union{MOI.GreaterThan{Float64},MOI.LessThan{Float64},MOI.EqualTo{Float64}}

mutable struct _VariableInfo
    index::MOI.VariableIndex
    column::Int
    bound::_BoundType
    # Both fields below are cached values to avoid triggering a model_update!
    # if the variable bounds are queried. They are NaN only if `bound` is
    # _NONE. _EQUAL_TO sets both of them. See also `lower_bound_if_soc`.
    lower_bound_if_bounded::Float64
    upper_bound_if_bounded::Float64
    type::Char
    start::Union{Float64,Nothing}
    name::String
    # Storage for constraint names associated with variables because Gurobi
    # can only store names for variables and proper constraints.
    # We can perform an optimization and only store three strings for the
    # constraint names because, at most, there can be three VariableIndex
    # constraints, e.g., LessThan, GreaterThan, and Integer.
    lessthan_name::String
    greaterthan_interval_or_equalto_name::String
    type_constraint_name::String
    # Storage for the lower bound if the variable is the `t` variable in a
    # second order cone. Theoretically, if both `lower_bound_if_bounded` and
    # `lower_bound_if_soc` are non-NaN, then they have the same value,
    # but you can also: (1) just have SOC constraints; (2) just have bounds;
    # (3) have a bound and a SOC constraint that does not need to set
    # `lower_bound_if_soc` (in all such cases just one of them is NaN).
    lower_bound_if_soc::Float64
    num_soc_constraints::Int
    function _VariableInfo(index::MOI.VariableIndex, column::Int)
        return new(
            index,
            column,
            _NONE,
            NaN,
            NaN,
            GRB_CONTINUOUS,
            nothing,
            "",
            "",
            "",
            "",
            NaN,
            0,
        )
    end
end

mutable struct _ConstraintInfo
    row::Int
    set::MOI.AbstractSet
    # Storage for constraint names. Where possible, these are also stored in
    # the Gurobi model.
    name::String
    _ConstraintInfo(row::Int, set) = new(row, set, "")
end

mutable struct Env
    ptr_env::Ptr{Cvoid}
    # These fields keep track of how many models the `Env` is used for to help
    # with finalizing. If you finalize an Env first, then the model, Gurobi will
    # throw an error.
    finalize_called::Bool
    attached_models::Int

    function Env(; output_flag::Int = 1)
        a = Ref{Ptr{Cvoid}}()
        ret = GRBemptyenv(a)
        env = new(a[], false, 0)
        _check_ret(env, ret)
        ret = GRBsetintparam(env.ptr_env, GRB_INT_PAR_OUTPUTFLAG, output_flag)
        _check_ret(env, ret)
        ret = GRBstartenv(env.ptr_env) 
        finalizer(env) do e
            e.finalize_called = true
            if e.attached_models == 0
                # Only finalize the model if there are no models using it.
                GRBfreeenv(e.ptr_env)
                e.ptr_env = C_NULL
            end
        end
        # Even if the loadenv fails, the pointer is still valid.
        _check_ret(env, ret)
        return env
    end
end

Base.cconvert(::Type{Ptr{Cvoid}}, x::Env) = x
Base.unsafe_convert(::Type{Ptr{Cvoid}}, env::Env) = env.ptr_env::Ptr{Cvoid}

const _HASH = CleverDicts.key_to_index
const _INVERSE_HASH = x -> CleverDicts.index_to_key(MOI.VariableIndex, x)

mutable struct Optimizer <: MOI.AbstractOptimizer
    # The low-level Gurobi model.
    inner::Ptr{Cvoid}
    # The Gurobi environment. If `nothing`, a new environment will be created
    # on `MOI.empty!`.
    env::Union{Nothing,Env}
    # The current user-provided parameters for the model.
    params::Dict{String,Any}

    # The next field is used to cleverly manage calls to `update_model!`.
    # `needs_update` is used to record whether an update should be called
    # before accessing a model attribute (such as the value of a RHS term).
    needs_update::Bool

    # A flag to keep track of MOI.Silent, which over-rides the OutputFlag
    # parameter.
    silent::Bool

    # An enum to remember what objective is currently stored in the model.
    objective_type::_ObjectiveType

    # track whether objective function is set and the state of objective sense
    is_objective_set::Bool
    objective_sense::Union{Nothing,MOI.OptimizationSense}

    # A mapping from the MOI.VariableIndex to the Gurobi column. _VariableInfo
    # also stores some additional fields like what bounds have been added, the
    # variable type, and the names of VariableIndex-in-Set constraints.
    variable_info::CleverDicts.CleverDict{
        MOI.VariableIndex,
        _VariableInfo,
        typeof(_HASH),
        typeof(_INVERSE_HASH),
    }

    # If you add variables to a model that had variables deleted AND has
    # not called `update_model!` since the deletion, then the newly created
    # variables may have attributes set, but their column index before the
    # call to `update_model!` is different than after the `update_model!`.
    # Before the `update_model!` their column is the same as if no variables
    # were deleted, after the `update_model!` the columns indexes are
    # shifted (by being being subtracted by the number of variables deleted
    # with column indexes smaller than them). To control this the two
    # fields below are used:
    # `next_column`: The column index of the next variable/column added. It is
    # updated when variables are added, and when the `_update_if_necessary!` is
    # called AND `columns_deleted_since_last_update` is not empty.
    # `columns_deleted_since_last_update`: Stores the column indexes of all
    # columns that were deleted since the last call to `_update_if_necessary!`,
    # after such call the vector is emptied.
    next_column::Int
    columns_deleted_since_last_update::Vector{Int}

    # An index that is incremented for each new constraint (regardless of type).
    # We can check if a constraint is valid by checking if it is in the correct
    # xxx_constraint_info. We should _not_ reset this to zero, since then new
    # constraints cannot be distinguished from previously created ones.
    last_constraint_index::Int
    # ScalarAffineFunction{Float64}-in-Set storage.
    affine_constraint_info::Dict{Int,_ConstraintInfo}
    # ScalarQuadraticFunction{Float64}-in-Set storage.
    quadratic_constraint_info::Dict{Int,_ConstraintInfo}
    # VectorOfVariables-in-Set storage.
    sos_constraint_info::Dict{Int,_ConstraintInfo}
    # VectorAffineFunction-in-Set storage.
    indicator_constraint_info::Dict{Int,_ConstraintInfo}
    # Note: we do not have a singlevariable_constraint_info dictionary. Instead,
    # data associated with these constraints are stored in the _VariableInfo
    # objects.

    # Mappings from variable and constraint names to their indices. These are
    # lazily built on-demand, so most of the time, they are `nothing`.
    name_to_variable::Union{
        Nothing,
        Dict{String,Union{Nothing,MOI.VariableIndex}},
    }
    name_to_constraint_index::Union{
        Nothing,
        Dict{String,Union{Nothing,MOI.ConstraintIndex}},
    }

    # Gurobi does not have a configurable memory limit (different of time),
    # but it does detect when it needs more memory than it is available,
    # and it stops the optimization returning a specific error code.
    # This is a different mechanism than Gurobi "Status" (that is used for
    # reporting why an optimization finished) and, in fact, may be triggered in
    # other cases than optimization (for example, when assembling the model).
    # For convenience, and homogeinity with other solvers, we save the code
    # returned by `GRBoptimize` in `ret_GRBoptimize`, and do not throw
    # an exception case it should be interpreted as a termination status.
    # Then, when/if the termination status is queried, we may override the
    # result taking into account the `ret_GRBoptimize` field.
    ret_GRBoptimize::Cint

    # These two flags allow us to distinguish between FEASIBLE_POINT and
    # INFEASIBILITY_CERTIFICATE when querying VariablePrimal and ConstraintDual.
    has_unbounded_ray::Bool
    has_infeasibility_cert::Bool

    # Callback fields.
    enable_interrupts::Bool
    callback_variable_primal::Vector{Float64}
    has_generic_callback::Bool
    callback_state::_CallbackState
    lazy_callback::Union{Nothing,Function}
    user_cut_callback::Union{Nothing,Function}
    heuristic_callback::Union{Nothing,Function}
    generic_callback::Any

    conflict::Cint

    """
        Optimizer(
            env::Union{Nothing,Env} = nothing;
            enable_interrupts::Bool = true,
        )

    Create a new Optimizer object.

    You can share Gurobi `Env`s between models by passing an instance of `Env`
    as the first argument.

    In order to enable interrupts via `CTRL+C`, a no-op callback is added to the
    model by default. In most cases, this has negligible effect on solution
    times. However, you can disable it (at the cost of not being able to
    interrupt a solve) by passing `enable_interrupts = false`.

    Set optimizer attributes using `MOI.RawOptimizerAttribute` or
    `JuMP.set_optimizer_atttribute`.

    ## Example

        using JuMP, Gurobi
        const env = Gurobi.Env()
        model = JuMP.Model(() -> Gurobi.Optimizer(env; enable_interrupts=false))
        set_optimizer_attribute(model, "OutputFlag", 0)
    """
    function Optimizer(
        env::Union{Nothing,Env} = nothing;
        enable_interrupts::Bool = true,
        kwargs...,
    )
        model = new()
        model.inner = C_NULL
        model.env = env === nothing ? Env() : env
        model.enable_interrupts = enable_interrupts
        model.params = Dict{String,Any}()
        if length(kwargs) > 0
            @warn("""Passing optimizer attributes as keyword arguments to
            Gurobi.Optimizer is deprecated. Use
                MOI.set(model, MOI.RawOptimizerAttribute("key"), value)
            or
                JuMP.set_optimizer_attribute(model, "key", value)
            instead.
            """)
        end
        for (name, value) in kwargs
            model.params[string(name)] = value
        end
        model.silent = false
        model.variable_info =
            CleverDicts.CleverDict{MOI.VariableIndex,_VariableInfo}(
                _HASH,
                _INVERSE_HASH,
            )
        model.next_column = 1
        model.last_constraint_index = 1
        model.columns_deleted_since_last_update = Int[]
        model.affine_constraint_info = Dict{Int,_ConstraintInfo}()
        model.quadratic_constraint_info = Dict{Int,_ConstraintInfo}()
        model.sos_constraint_info = Dict{Int,_ConstraintInfo}()
        model.indicator_constraint_info = Dict{Int,_ConstraintInfo}()

        model.callback_variable_primal = Float64[]
        MOI.empty!(model)
        finalizer(model) do m
            ret = GRBfreemodel(m.inner)
            _check_ret(m, ret)
            m.env.attached_models -= 1
            if env === nothing
                @assert m.env.attached_models == 0
                # We created this environment. Finalize it now.
                finalize(m.env)
            elseif m.env.finalize_called && m.env.attached_models == 0
                # We delayed finalizing `m.env` earlier because there were still
                # models attached. Finalize it now.
                GRBfreeenv(m.env.ptr_env)
                m.env.ptr_env = C_NULL
            end
        end
        return model
    end
end
Base.cconvert(::Type{Ptr{Cvoid}}, x::Optimizer) = x
function Base.unsafe_convert(::Type{Ptr{Cvoid}}, model::Optimizer)
    return model.inner::Ptr{Cvoid}
end

function _check_ret(model::Optimizer, ret::Cint)
    if ret != 0
        msg = unsafe_string(GRBgetmerrormsg(model))
        throw(ErrorException("Gurobi Error $(ret): $(msg)"))
    end
    return
end

# If you add a new error code that, when returned by GRBoptimize,
# should be treated as a TerminationStatus by MOI, to the global `Dict`
# below, then the rest of the code should pick up on this seamlessly.
const _ERROR_TO_STATUS = Dict{Cint,Tuple{MOI.TerminationStatusCode,String}}([
    # Code => (TerminationStatus, RawStatusString)
    GRB_ERROR_OUT_OF_MEMORY =>
        (MOI.MEMORY_LIMIT, "Available memory was exhausted."),
])

# Same as _check_ret, but deals with the `model.ret_GRBoptimize` machinery.
function _check_ret_GRBoptimize(model)
    if !haskey(_ERROR_TO_STATUS, model.ret_GRBoptimize)
        _check_ret(model, model.ret_GRBoptimize)
    end
    return
end

function _check_ret(env, ret::Cint)
    if ret != 0
        msg = unsafe_string(GRBgeterrormsg(env))
        throw(ErrorException("Gurobi Error $(ret): $(msg)"))
    end
    return
end

function Base.show(io::IO, model::Optimizer)
    if model.inner == C_NULL
        println(io, "Gurobi Model: NULL")
        return
    end
    p = Ref{Cint}()
    GRBgetintattr(model, "ModelSense", p)
    println(io, "    sense  : $(p[] > 0 ? :minimize : :maximize)")
    GRBgetintattr(model, "NumVars", p)
    println(io, "    number of variables             = $(p[])")
    GRBgetintattr(model, "NumConstrs", p)
    println(io, "    number of linear constraints    = $(p[])")
    GRBgetintattr(model, "NumQConstrs", p)
    println(io, "    number of quadratic constraints = $(p[])")
    GRBgetintattr(model, "NumSOS", p)
    println(io, "    number of sos constraints       = $(p[])")
    GRBgetintattr(model, "NumNZs", p)
    println(io, "    number of non-zero coeffs       = $(p[])")
    GRBgetintattr(model, "NumQNZs", p)
    println(io, "    number of non-zero qp objective terms  = $(p[])")
    GRBgetintattr(model, "NumQCNZs", p)
    println(io, "    number of non-zero qp constraint terms = $(p[])")
    return
end

function MOI.empty!(model::Optimizer)
    # Free the current model, if it exists.
    if model.inner != C_NULL
        ret = GRBfreemodel(model.inner)
        _check_ret(model, ret)
        model.env.attached_models -= 1
    end
    # Then create a new one
    a = Ref{Ptr{Cvoid}}()
    ret =
        GRBnewmodel(model.env, a, "", 0, C_NULL, C_NULL, C_NULL, C_NULL, C_NULL)
    model.inner = a[]
    model.env.attached_models += 1
    _check_ret(model, ret)
    # Reset the parameters in this new environment
    if model.silent
        MOI.set(model, MOI.Silent(), true)
    end
    for (name, value) in model.params
        MOI.set(model, MOI.RawOptimizerAttribute(name), value)
    end
    model.needs_update = false
    model.objective_type = _SCALAR_AFFINE
    model.is_objective_set = false
    model.objective_sense = nothing
    empty!(model.variable_info)
    model.next_column = 1
    empty!(model.columns_deleted_since_last_update)
    empty!(model.affine_constraint_info)
    empty!(model.quadratic_constraint_info)
    empty!(model.sos_constraint_info)
    empty!(model.indicator_constraint_info)
    model.name_to_variable = nothing
    model.name_to_constraint_index = nothing
    model.ret_GRBoptimize = Cint(0)
    model.has_unbounded_ray = false
    model.has_infeasibility_cert = false
    empty!(model.callback_variable_primal)
    model.callback_state = _CB_NONE
    model.has_generic_callback = false
    model.lazy_callback = nothing
    model.user_cut_callback = nothing
    model.heuristic_callback = nothing
    model.generic_callback = nothing
    model.conflict = Cint(-1)
    return
end

function MOI.is_empty(model::Optimizer)
    model.needs_update && return false
    model.objective_type != _SCALAR_AFFINE && return false
    model.is_objective_set == true && return false
    model.objective_sense !== nothing && return false
    !isempty(model.variable_info) && return false
    !isone(model.next_column) && return false
    !isempty(model.columns_deleted_since_last_update) && return false
    !isempty(model.affine_constraint_info) && return false
    !isempty(model.quadratic_constraint_info) && return false
    !isempty(model.sos_constraint_info) && return false
    model.name_to_variable !== nothing && return false
    model.name_to_constraint_index !== nothing && return false
    !iszero(model.ret_GRBoptimize) && return false
    model.has_unbounded_ray && return false
    model.has_infeasibility_cert && return false
    !isempty(model.callback_variable_primal) && return false
    model.callback_state != _CB_NONE && return false
    model.has_generic_callback && return false
    model.lazy_callback !== nothing && return false
    model.user_cut_callback !== nothing && return false
    model.heuristic_callback !== nothing && return false
    return true
end

"""
    _require_update(model::Optimizer)

Sets the `model.needs_update` flag. Call this at the end of any mutating method.
"""
function _require_update(model::Optimizer)
    model.needs_update = true
    return
end

"""
    _update_if_necessary(model::Optimizer)

Calls `update_model!`, but only if the `model.needs_update` flag is set.
"""
function _update_if_necessary(model::Optimizer)
    if model.needs_update
        sort!(model.columns_deleted_since_last_update)
        for var_info in values(model.variable_info)
            # The trick here is: searchsortedlast returns, in O(log n), the
            # last index with a column smaller than var_info.column, over
            # columns_deleted_since_last_update this is the same as the number
            # of columns deleted before it, and how much its value need to be
            # shifted.
            var_info.column -= searchsortedlast(
                model.columns_deleted_since_last_update,
                var_info.column,
            )
        end
        model.next_column -= length(model.columns_deleted_since_last_update)
        empty!(model.columns_deleted_since_last_update)
        ret = GRBupdatemodel(model)
        _check_ret(model, ret)
        model.needs_update = false
    else
        @assert isempty(model.columns_deleted_since_last_update)
    end
    return
end

MOI.get(::Optimizer, ::MOI.SolverName) = "Gurobi"

MOI.get(::Optimizer, ::MOI.SolverVersion) = string(_GUROBI_VERSION)

function MOI.supports(
    ::Optimizer,
    ::MOI.ObjectiveFunction{F},
) where {
    F<:Union{
        MOI.VariableIndex,
        MOI.ScalarAffineFunction{Float64},
        MOI.ScalarQuadraticFunction{Float64},
    },
}
    return true
end

function MOI.supports_constraint(
    ::Optimizer,
    ::Type{MOI.VariableIndex},
    ::Type{F},
) where {
    F<:Union{
        MOI.EqualTo{Float64},
        MOI.LessThan{Float64},
        MOI.GreaterThan{Float64},
        MOI.Interval{Float64},
        MOI.ZeroOne,
        MOI.Integer,
        MOI.Semicontinuous{Float64},
        MOI.Semiinteger{Float64},
    },
}
    return true
end

function MOI.supports_constraint(
    ::Optimizer,
    ::Type{MOI.VectorOfVariables},
    ::Type{F},
) where {F<:Union{MOI.SOS1{Float64},MOI.SOS2{Float64},MOI.SecondOrderCone}}
    return true
end

# We choose _not_ to support ScalarAffineFunction-in-Interval and
# ScalarQuadraticFunction-in-Interval because Gurobi introduces some slack
# variables that makes it hard to keep track of the column indices.

function MOI.supports_constraint(
    ::Optimizer,
    ::Type{MOI.ScalarAffineFunction{Float64}},
    ::Type{F},
) where {
    F<:Union{
        MOI.EqualTo{Float64},
        MOI.LessThan{Float64},
        MOI.GreaterThan{Float64},
    },
}
    return true
end

function MOI.supports_constraint(
    ::Optimizer,
    ::Type{MOI.ScalarQuadraticFunction{Float64}},
    ::Type{F},
) where {
    F<:Union{
        MOI.EqualTo{Float64},
        MOI.LessThan{Float64},
        MOI.GreaterThan{Float64},
    },
}
    return true
end

MOI.supports(::Optimizer, ::MOI.VariableName, ::Type{MOI.VariableIndex}) = true
function MOI.supports(
    ::Optimizer,
    ::MOI.ConstraintName,
    ::Type{<:MOI.ConstraintIndex},
)
    return true
end

MOI.supports(::Optimizer, ::MOI.Name) = true
MOI.supports(::Optimizer, ::MOI.Silent) = true
MOI.supports(::Optimizer, ::MOI.NumberOfThreads) = true
MOI.supports(::Optimizer, ::MOI.TimeLimitSec) = true
MOI.supports(::Optimizer, ::MOI.ObjectiveSense) = true
MOI.supports(::Optimizer, ::MOI.RawOptimizerAttribute) = true
MOI.supports(::Optimizer, ::MOI.ConstraintPrimalStart) = false
MOI.supports(::Optimizer, ::MOI.ConstraintDualStart) = false

function MOI.set(model::Optimizer, raw::MOI.RawOptimizerAttribute, value)
    env = GRBgetenv(model)
    param = raw.name
    model.params[param] = value
    param_type = GRBgetparamtype(env, param)
    ret = if param_type == -1
        throw(MOI.UnsupportedAttribute(MOI.RawOptimizerAttribute(param)))
    elseif param_type == 1
        GRBsetintparam(env, param, value)
    elseif param_type == 2
        GRBsetdblparam(env, param, value)
    else
        @assert param_type == 3
        GRBsetstrparam(env, param, value)
    end
    return _check_ret(env, ret)
end

function MOI.get(model::Optimizer, raw::MOI.RawOptimizerAttribute)
    env = GRBgetenv(model)
    param = raw.name
    param_type = GRBgetparamtype(env, param)
    if param_type == -1
        throw(MOI.UnsupportedAttribute(MOI.RawOptimizerAttribute(param)))
    elseif param_type == 1
        a = Ref{Cint}()
        ret = GRBgetintparam(env, param, a)
        _check_ret(env, ret)
        return a[]
    elseif param_type == 2
        a = Ref{Cdouble}()
        ret = GRBgetdblparam(env, param, a)
        _check_ret(env, ret)
        return a[]
    else
        @assert param_type == 3
        valueP = Ref{Ptr{Cchar}}()
        ret = GRBgetstrparam(env, param, valueP)
        _check_ret(env, ret)
        return unsafe_string(valueP[])
    end
end

function MOI.set(model::Optimizer, ::MOI.TimeLimitSec, limit::Real)
    MOI.set(model, MOI.RawOptimizerAttribute("TimeLimit"), limit)
    return
end

function MOI.get(model::Optimizer, ::MOI.TimeLimitSec)
    return MOI.get(model, MOI.RawOptimizerAttribute("TimeLimit"))
end

MOI.supports_incremental_interface(::Optimizer) = true

function MOI.copy_to(dest::Optimizer, src::MOI.ModelLike)
    return MOI.Utilities.default_copy_to(dest, src)
end

function MOI.get(model::Optimizer, ::MOI.ListOfVariableAttributesSet)
    ret = MOI.AbstractVariableAttribute[]
    found_name, found_start = false, false
    for info in values(model.variable_info)
        if !found_name && !isempty(info.name)
            push!(ret, MOI.VariableName())
            found_name = true
        end
        if !found_start && info.start !== nothing
            push!(ret, MOI.VariablePrimalStart())
            found_start = true
        end
        if found_start && found_name
            return ret
        end
    end
    return ret
end

function MOI.get(model::Optimizer, ::MOI.ListOfModelAttributesSet)
    if MOI.is_empty(model)
        return Any[]
    end
    attributes = Any[]
    if model.objective_sense !== nothing
        push!(attributes, MOI.ObjectiveSense())
    end
    if model.is_objective_set
        F = MOI.get(model, MOI.ObjectiveFunctionType())
        push!(attributes, MOI.ObjectiveFunction{F}())
    end
    if MOI.get(model, MOI.Name()) != ""
        push!(attributes, MOI.Name())
    end
    return attributes
end

function MOI.get(
    model::Optimizer,
    ::MOI.ListOfConstraintAttributesSet{F,S},
) where {S,F}
    ret = MOI.AbstractConstraintAttribute[]
    constraint_indices = MOI.get(model, MOI.ListOfConstraintIndices{F,S}())
    found_name = any(
        !isempty(MOI.get(model, MOI.ConstraintName(), index)) for
        index in constraint_indices
    )
    if found_name
        push!(ret, MOI.ConstraintName())
    end
    return ret
end

function _indices_and_coefficients(
    indices::AbstractVector{Cint},
    coefficients::AbstractVector{Float64},
    model::Optimizer,
    f::MOI.ScalarAffineFunction{Float64},
)
    i = 1
    for term in f.terms
        indices[i] = Cint(column(model, term.variable) - 1)
        coefficients[i] = term.coefficient
        i += 1
    end
    return indices, coefficients
end

function _indices_and_coefficients(
    model::Optimizer,
    f::MOI.ScalarAffineFunction{Float64},
)
    f_canon = if MOI.Utilities.is_canonical(f)
        f
    else
        MOI.Utilities.canonical(f)
    end
    nnz = length(f_canon.terms)
    indices = Vector{Cint}(undef, nnz)
    coefficients = Vector{Float64}(undef, nnz)
    _indices_and_coefficients(indices, coefficients, model, f_canon)
    return indices, coefficients
end

function _indices_and_coefficients(
    I::AbstractVector{Cint},
    J::AbstractVector{Cint},
    V::AbstractVector{Float64},
    indices::AbstractVector{Cint},
    coefficients::AbstractVector{Float64},
    model::Optimizer,
    f::MOI.ScalarQuadraticFunction,
)
    for (i, term) in enumerate(f.quadratic_terms)
        I[i] = Cint(column(model, term.variable_1) - 1)
        J[i] = Cint(column(model, term.variable_2) - 1)
        V[i] = term.coefficient
        # Gurobi returns a list of terms. MOI requires 0.5 x' Q x. So, to get
        # from
        #   Gurobi -> MOI => multiply diagonals by 2.0
        #   MOI -> Gurobi => multiply diagonals by 0.5
        # Example: 2x^2 + x*y + y^2
        #   |x y| * |a b| * |x| = |ax+by bx+cy| * |x| = 0.5ax^2 + bxy + 0.5cy^2
        #           |b c|   |y|                   |y|
        #   Gurobi needs: (I, J, V) = ([0, 0, 1], [0, 1, 1], [2, 1, 1])
        #   MOI needs:
        #     [SQT(4.0, x, x), SQT(1.0, x, y), SQT(2.0, y, y)]
        if I[i] == J[i]
            V[i] *= 0.5
        end
    end
    for (i, term) in enumerate(f.affine_terms)
        indices[i] = Cint(column(model, term.variable) - 1)
        coefficients[i] = term.coefficient
    end
    return
end

function _indices_and_coefficients(
    model::Optimizer,
    f::MOI.ScalarQuadraticFunction,
)
    f_canon = if MOI.Utilities.is_canonical(f)
        f
    else
        MOI.Utilities.canonical(f)
    end
    nnz_quadratic = length(f_canon.quadratic_terms)
    nnz_affine = length(f_canon.affine_terms)
    I = Vector{Cint}(undef, nnz_quadratic)
    J = Vector{Cint}(undef, nnz_quadratic)
    V = Vector{Float64}(undef, nnz_quadratic)
    indices = Vector{Cint}(undef, nnz_affine)
    coefficients = Vector{Float64}(undef, nnz_affine)
    _indices_and_coefficients(I, J, V, indices, coefficients, model, f_canon)
    return indices, coefficients, I, J, V
end

_sense_and_rhs(s::MOI.LessThan{Float64}) = (GRB_LESS_EQUAL, s.upper)
_sense_and_rhs(s::MOI.GreaterThan{Float64}) = (GRB_GREATER_EQUAL, s.lower)
_sense_and_rhs(s::MOI.EqualTo{Float64}) = (GRB_EQUAL, s.value)

###
### Variables
###

# Short-cuts to return the _VariableInfo associated with an index.
function _info(model::Optimizer, key::MOI.VariableIndex)
    if haskey(model.variable_info, key)
        return model.variable_info[key]
    end
    return throw(MOI.InvalidIndex(key))
end

"""
    column(model::Optimizer, x::MOI.VariableIndex)

Return the 1-indexed column associated with `x`.

The C API requires 0-indexed columns.
"""
function column(model::Optimizer, x::MOI.VariableIndex)
    return _info(model, x).column
end

function _get_next_column(model::Optimizer)
    model.next_column += 1
    return model.next_column - 1
end

function MOI.add_variable(model::Optimizer)
    # Initialize `_VariableInfo` with a dummy `VariableIndex` and a column,
    # because we need `add_item` to tell us what the `VariableIndex` is.
    index = CleverDicts.add_item(
        model.variable_info,
        _VariableInfo(MOI.VariableIndex(0), 0),
    )
    info = _info(model, index)
    # Now, set `.index` and `.column`.
    info.index = index
    info.column = _get_next_column(model)
    ret =
        GRBaddvar(model, 0, C_NULL, C_NULL, 0.0, -Inf, Inf, GRB_CONTINUOUS, "")
    _check_ret(model, ret)
    _require_update(model)
    return index
end

function MOI.add_variables(model::Optimizer, N::Int)
    ret = GRBaddvars(
        model,
        N,
        0,
        C_NULL,
        C_NULL,
        C_NULL,
        C_NULL,
        fill(-Inf, N),
        C_NULL,
        C_NULL,
        C_NULL,
    )
    _check_ret(model, ret)
    indices = Vector{MOI.VariableIndex}(undef, N)
    for i in 1:N
        # Initialize `_VariableInfo` with a dummy `VariableIndex` and a column,
        # because we need `add_item` to tell us what the `VariableIndex` is.
        index = CleverDicts.add_item(
            model.variable_info,
            _VariableInfo(MOI.VariableIndex(0), 0),
        )
        info = _info(model, index)
        # Now, set `.index` and `.column`.
        info.index = index
        info.column = _get_next_column(model)
        indices[i] = index
    end
    _require_update(model)
    return indices
end

# We implement a specialized version here to avoid calling into Gurobi twice.
# Using the standard implementation, we would first create a variable in Gurobi
# with GRBaddvar that has bounds of (-Inf,+Inf), and then immediately after
# reset those bounds using the attributes interface. Instead, we just pass the
# desired bounds directly to GRBaddvar.
function MOI.add_constrained_variable(
    model::Optimizer,
    set::S,
)::Tuple{
    MOI.VariableIndex,
    MOI.ConstraintIndex{MOI.VariableIndex,S},
} where {S<:_SCALAR_SETS}
    vi = CleverDicts.add_item(
        model.variable_info,
        _VariableInfo(MOI.VariableIndex(0), 0),
    )
    info = _info(model, vi)
    # Now, set `.index` and `.column`.
    info.index = vi
    info.column = _get_next_column(model)
    lb = -Inf
    ub = Inf
    if S <: MOI.LessThan{Float64}
        ub = set.upper
        info.upper_bound_if_bounded = ub
        info.bound = _LESS_THAN
    elseif S <: MOI.GreaterThan{Float64}
        lb = set.lower
        info.lower_bound_if_bounded = lb
        info.bound = _GREATER_THAN
    elseif S <: MOI.EqualTo{Float64}
        lb = set.value
        ub = set.value
        info.lower_bound_if_bounded = lb
        info.upper_bound_if_bounded = ub
        info.bound = _EQUAL_TO
    else
        @assert S <: MOI.Interval{Float64}
        lb = set.lower
        ub = set.upper
        info.lower_bound_if_bounded = lb
        info.upper_bound_if_bounded = ub
        info.bound = _INTERVAL
    end
    ret = GRBaddvar(model, 0, C_NULL, C_NULL, 0.0, lb, ub, GRB_CONTINUOUS, "")
    _check_ret(model, ret)
    _require_update(model)
    ci = MOI.ConstraintIndex{MOI.VariableIndex,typeof(set)}(vi.value)
    return vi, ci
end

function MOI.is_valid(model::Optimizer, v::MOI.VariableIndex)
    return haskey(model.variable_info, v)
end

function MOI.delete(model::Optimizer, indices::Vector{<:MOI.VariableIndex})
    #_update_if_necessary(model)
    info = [_info(model, var_idx) for var_idx in indices]
    soc_idx = findfirst(e -> e.num_soc_constraints > 0, info)
    soc_idx !== nothing && throw(MOI.DeleteNotAllowed(indices[soc_idx]))
    del_cols = collect(Cint(i.column - 1) for i in info)
    ret = GRBdelvars(model, length(del_cols), del_cols)
    _check_ret(model, ret)
    for var_idx in indices
        delete!(model.variable_info, var_idx)
    end
    append!(model.columns_deleted_since_last_update, del_cols .+ 1)
    model.name_to_variable = nothing
    # We throw away name_to_constraint_index so we will rebuild VariableIndex
    # constraint names without v.
    model.name_to_constraint_index = nothing
    _require_update(model)
    return
end

function MOI.delete(model::Optimizer, v::MOI.VariableIndex)
    #_update_if_necessary(model)
    info = _info(model, v)
    if info.num_soc_constraints > 0
        throw(MOI.DeleteNotAllowed(v))
    end
    push!(model.columns_deleted_since_last_update, info.column)
    ret = GRBdelvars(model, 1, Ref{Cint}(info.column - 1))
    _check_ret(model, ret)
    delete!(model.variable_info, v)
    model.name_to_variable = nothing
    # We throw away name_to_constraint_index so we will rebuild VariableIndex
    # constraint names without v.
    model.name_to_constraint_index = nothing
    _require_update(model)
    return
end

function MOI.get(model::Optimizer, ::Type{MOI.VariableIndex}, name::String)
    if model.name_to_variable === nothing
        _rebuild_name_to_variable(model)
    end
    if haskey(model.name_to_variable, name)
        variable = model.name_to_variable[name]
        if variable === nothing
            error("Duplicate variable name detected: $(name)")
        end
        return variable
    end
    return nothing
end

function _rebuild_name_to_variable(model::Optimizer)
    model.name_to_variable = Dict{String,Union{Nothing,MOI.VariableIndex}}()
    for (index, info) in model.variable_info
        if isempty(info.name)
            continue
        end
        if haskey(model.name_to_variable, info.name)
            model.name_to_variable[info.name] = nothing
        else
            model.name_to_variable[info.name] = index
        end
    end
    return
end

function MOI.get(model::Optimizer, ::MOI.VariableName, v::MOI.VariableIndex)
    return _info(model, v).name
end

function MOI.set(
    model::Optimizer,
    ::MOI.VariableName,
    v::MOI.VariableIndex,
    name::String,
)
    info = _info(model, v)
    info.name = name
    ret = GRBsetstrattrelement(model, "VarName", Cint(info.column) - 1, name)
    _check_ret(model, ret)
    _require_update(model)
    model.name_to_variable = nothing
    return
end

###
### Objectives
###

function _zero_objective(model::Optimizer)
    num_vars = length(model.variable_info)
    obj = zeros(Float64, num_vars)
    _update_if_necessary(model)
    ret = GRBdelq(model)
    _check_ret(model, ret)
    ret = GRBsetdblattrarray(model, "Obj", 0, num_vars, obj)
    _check_ret(model, ret)
    ret = GRBsetdblattr(model, "ObjCon", 0.0)
    _check_ret(model, ret)
    return _require_update(model)
end

function MOI.set(
    model::Optimizer,
    ::MOI.ObjectiveSense,
    sense::MOI.OptimizationSense,
)
    if sense == MOI.MIN_SENSE
        ret = GRBsetintattr(model, "ModelSense", 1)
        _check_ret(model, ret)
    elseif sense == MOI.MAX_SENSE
        ret = GRBsetintattr(model, "ModelSense", -1)
        _check_ret(model, ret)
    else
        @assert sense == MOI.FEASIBILITY_SENSE
        _zero_objective(model)
        ret = GRBsetintattr(model, "ModelSense", 1)
        _check_ret(model, ret)
    end
    model.objective_sense = sense
    _require_update(model)
    return
end

function MOI.get(model::Optimizer, ::MOI.ObjectiveSense)
    if model.objective_sense !== nothing
        return model.objective_sense
    else
        return MOI.FEASIBILITY_SENSE
    end
end

function MOI.set(
    model::Optimizer,
    ::MOI.ObjectiveFunction{F},
    f::F,
) where {F<:MOI.VariableIndex}
    MOI.set(
        model,
        MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(),
        convert(MOI.ScalarAffineFunction{Float64}, f),
    )
    model.objective_type = _SINGLE_VARIABLE
    model.is_objective_set = true
    return
end

function MOI.get(model::Optimizer, ::MOI.ObjectiveFunction{MOI.VariableIndex})
    obj = MOI.get(
        model,
        MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(),
    )
    return convert(MOI.VariableIndex, obj)
end

function MOI.set(
    model::Optimizer,
    ::MOI.ObjectiveFunction{F},
    f::F,
) where {F<:MOI.ScalarAffineFunction{Float64}}
    if model.objective_type == _SCALAR_QUADRATIC
        # We need to zero out the existing quadratic objective.
        ret = GRBdelq(model)
        _check_ret(model, ret)
    end
    num_vars = length(model.variable_info)
    obj = zeros(Float64, num_vars)
    for term in f.terms
        obj[column(model, term.variable)] += term.coefficient
    end
    # NOTE: variables added may be referred before a `_update_if_necessary`
    # what is the problem we try to prevent below?
    # This update is needed because we might have added some variables.
    _update_if_necessary(model)
    ret = GRBsetdblattrarray(model, "Obj", 0, num_vars, obj)
    _check_ret(model, ret)
    ret = GRBsetdblattr(model, "ObjCon", f.constant)
    _check_ret(model, ret)
    _require_update(model)
    model.objective_type = _SCALAR_AFFINE
    model.is_objective_set = true
    return
end

function MOI.get(
    model::Optimizer,
    ::MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}},
)
    if model.objective_type == _SCALAR_QUADRATIC
        error(
            "Unable to get objective function. Currently: " *
            "$(model.objective_type).",
        )
    end
    _update_if_necessary(model)
    dest = zeros(length(model.variable_info))
    ret = GRBgetdblattrarray(model, "Obj", 0, length(dest), dest)
    _check_ret(model, ret)
    terms = MOI.ScalarAffineTerm{Float64}[]
    for (index, info) in model.variable_info
        coefficient = dest[info.column]
        iszero(coefficient) && continue
        push!(terms, MOI.ScalarAffineTerm(coefficient, index))
    end
    constant = Ref{Cdouble}()
    ret = GRBgetdblattr(model, "ObjCon", constant)
    _check_ret(model, ret)
    return MOI.ScalarAffineFunction(terms, constant[])
end

function MOI.set(
    model::Optimizer,
    ::MOI.ObjectiveFunction{F},
    f::F,
) where {F<:MOI.ScalarQuadraticFunction{Float64}}
    affine_indices, affine_coefficients, I, J, V =
        _indices_and_coefficients(model, f)
    _update_if_necessary(model)
    # We need to zero out any existing linear objective.
    obj = zeros(length(model.variable_info))
    for (i, c) in zip(affine_indices, affine_coefficients)
        obj[i+1] = c
    end
    ret = GRBsetdblattrarray(model, "Obj", 0, length(obj), obj)
    _check_ret(model, ret)
    ret = GRBsetdblattr(model, "ObjCon", f.constant)
    _check_ret(model, ret)
    # We need to zero out the existing quadratic objective.
    ret = GRBdelq(model)
    _check_ret(model, ret)
    ret = GRBaddqpterms(model, length(I), I, J, V)
    _check_ret(model, ret)
    _require_update(model)
    model.objective_type = _SCALAR_QUADRATIC
    model.is_objective_set = true
    return
end

function MOI.get(
    model::Optimizer,
    ::MOI.ObjectiveFunction{MOI.ScalarQuadraticFunction{Float64}},
)
    _update_if_necessary(model)
    dest = zeros(length(model.variable_info))
    ret = GRBgetdblattrarray(model, "Obj", 0, length(dest), dest)
    _check_ret(model, ret)
    terms = MOI.ScalarAffineTerm{Float64}[]
    for (index, info) in model.variable_info
        coefficient = dest[info.column]
        iszero(coefficient) && continue
        push!(terms, MOI.ScalarAffineTerm(coefficient, index))
    end
    constant = Ref{Cdouble}()
    ret = GRBgetdblattr(model, "ObjCon", constant)
    _check_ret(model, ret)
    q_terms = MOI.ScalarQuadraticTerm{Float64}[]

    numqnzP = Ref{Cint}()
    ret = GRBgetintattr(model, "NumQNZs", numqnzP)
    _check_ret(model, ret)
    qrow = Array{Cint}(undef, numqnzP[])
    qcol = Array{Cint}(undef, numqnzP[])
    qval = Array{Float64}(undef, numqnzP[])
    nzout = Ref{Cint}()
    ret = GRBgetq(model, numqnzP, qrow, qcol, qval)
    _check_ret(model, ret)
    for (i, j, v) in zip(qrow, qcol, qval)
        if iszero(v)
            continue
        end
        # See note in `_indices_and_coefficients`.
        new_v = i == j ? 2v : v
        push!(
            q_terms,
            MOI.ScalarQuadraticTerm(
                new_v,
                model.variable_info[CleverDicts.LinearIndex(i + 1)].index,
                model.variable_info[CleverDicts.LinearIndex(j + 1)].index,
            ),
        )
    end
    return MOI.ScalarQuadraticFunction(q_terms, terms, constant[])
end

function MOI.modify(
    model::Optimizer,
    ::MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}},
    chg::MOI.ScalarConstantChange{Float64},
)
    ret = GRBsetdblattr(model, "ObjCon", chg.new_constant)
    _check_ret(model, ret)
    _require_update(model)
    return
end

##
##  VariableIndex-in-Set constraints.
##

function _info(
    model::Optimizer,
    c::MOI.ConstraintIndex{MOI.VariableIndex,<:Any},
)
    var_index = MOI.VariableIndex(c.value)
    if haskey(model.variable_info, var_index)
        return _info(model, var_index)
    end
    return throw(MOI.InvalidIndex(c))
end

"""
    column(model::Optimizer, c::MOI.ConstraintIndex{MOI.VariableIndex, <:Any})

Return the 1-indexed column associated with `c`.

The C API requires 0-indexed columns.
"""
function column(
    model::Optimizer,
    c::MOI.ConstraintIndex{MOI.VariableIndex,<:Any},
)
    return _info(model, c).column
end

function MOI.is_valid(
    model::Optimizer,
    c::MOI.ConstraintIndex{MOI.VariableIndex,MOI.LessThan{Float64}},
)
    if haskey(model.variable_info, MOI.VariableIndex(c.value))
        info = _info(model, c)
        return info.bound == _LESS_THAN || info.bound == _LESS_AND_GREATER_THAN
    end
    return false
end

function MOI.is_valid(
    model::Optimizer,
    c::MOI.ConstraintIndex{MOI.VariableIndex,MOI.GreaterThan{Float64}},
)
    if haskey(model.variable_info, MOI.VariableIndex(c.value))
        info = _info(model, c)
        return info.bound == _GREATER_THAN ||
               info.bound == _LESS_AND_GREATER_THAN
    end
    return false
end

function MOI.is_valid(
    model::Optimizer,
    c::MOI.ConstraintIndex{MOI.VariableIndex,MOI.Interval{Float64}},
)
    return haskey(model.variable_info, MOI.VariableIndex(c.value)) &&
           _info(model, c).bound == _INTERVAL
end

function MOI.is_valid(
    model::Optimizer,
    c::MOI.ConstraintIndex{MOI.VariableIndex,MOI.EqualTo{Float64}},
)
    return haskey(model.variable_info, MOI.VariableIndex(c.value)) &&
           _info(model, c).bound == _EQUAL_TO
end

function MOI.is_valid(
    model::Optimizer,
    c::MOI.ConstraintIndex{MOI.VariableIndex,MOI.ZeroOne},
)
    return haskey(model.variable_info, MOI.VariableIndex(c.value)) &&
           _info(model, c).type == GRB_BINARY
end

function MOI.is_valid(
    model::Optimizer,
    c::MOI.ConstraintIndex{MOI.VariableIndex,MOI.Integer},
)
    return haskey(model.variable_info, MOI.VariableIndex(c.value)) &&
           _info(model, c).type == GRB_INTEGER
end

function MOI.is_valid(
    model::Optimizer,
    c::MOI.ConstraintIndex{MOI.VariableIndex,MOI.Semicontinuous{Float64}},
)
    return haskey(model.variable_info, MOI.VariableIndex(c.value)) &&
           _info(model, c).type == GRB_SEMICONT
end

function MOI.is_valid(
    model::Optimizer,
    c::MOI.ConstraintIndex{MOI.VariableIndex,MOI.Semiinteger{Float64}},
)
    return haskey(model.variable_info, MOI.VariableIndex(c.value)) &&
           _info(model, c).type == GRB_SEMIINT
end

function MOI.get(
    model::Optimizer,
    ::MOI.ConstraintFunction,
    c::MOI.ConstraintIndex{MOI.VariableIndex,<:Any},
)
    MOI.throw_if_not_valid(model, c)
    return MOI.VariableIndex(c.value)
end

function MOI.set(
    model::Optimizer,
    ::MOI.ConstraintFunction,
    c::MOI.ConstraintIndex{MOI.VariableIndex,<:Any},
    ::MOI.VariableIndex,
)
    return throw(MOI.SettingVariableIndexNotAllowed())
end

_bounds(s::MOI.GreaterThan{Float64}) = (s.lower, nothing)
_bounds(s::MOI.LessThan{Float64}) = (nothing, s.upper)
_bounds(s::MOI.EqualTo{Float64}) = (s.value, s.value)
_bounds(s::MOI.Interval{Float64}) = (s.lower, s.upper)

function _throw_if_existing_lower(
    bound::_BoundType,
    var_type::Char,
    new_set::Type{<:MOI.AbstractSet},
    variable::MOI.VariableIndex,
)
    existing_set = if bound == _LESS_AND_GREATER_THAN || bound == _GREATER_THAN
        MOI.GreaterThan{Float64}
    elseif bound == _INTERVAL
        MOI.Interval{Float64}
    elseif bound == _EQUAL_TO
        MOI.EqualTo{Float64}
    elseif var_type == GRB_SEMIINT
        MOI.Semiinteger{Float64}
    elseif var_type == GRB_SEMICONT
        MOI.Semicontinuous{Float64}
    else
        nothing  # Also covers `_NONE` and `_LESS_THAN`.
    end
    if existing_set !== nothing
        throw(MOI.LowerBoundAlreadySet{existing_set,new_set}(variable))
    end
end

function _throw_if_existing_upper(
    bound::_BoundType,
    var_type::Char,
    new_set::Type{<:MOI.AbstractSet},
    variable::MOI.VariableIndex,
)
    existing_set = if bound == _LESS_AND_GREATER_THAN || bound == _LESS_THAN
        MOI.LessThan{Float64}
    elseif bound == _INTERVAL
        MOI.Interval{Float64}
    elseif bound == _EQUAL_TO
        MOI.EqualTo{Float64}
    elseif var_type == GRB_SEMIINT
        MOI.Semiinteger{Float64}
    elseif var_type == GRB_SEMICONT
        MOI.Semicontinuous{Float64}
    else
        nothing  # Also covers `_NONE` and `_GREATER_THAN`.
    end
    if existing_set !== nothing
        throw(MOI.UpperBoundAlreadySet{existing_set,new_set}(variable))
    end
end

function MOI.add_constraint(
    model::Optimizer,
    f::MOI.VariableIndex,
    s::S,
) where {S<:_SCALAR_SETS}
    info = _info(model, f)
    if S <: MOI.LessThan{Float64}
        _throw_if_existing_upper(info.bound, info.type, S, f)
        info.bound =
            info.bound == _GREATER_THAN ? _LESS_AND_GREATER_THAN : _LESS_THAN
        info.upper_bound_if_bounded = s.upper
    elseif S <: MOI.GreaterThan{Float64}
        _throw_if_existing_lower(info.bound, info.type, S, f)
        info.bound =
            info.bound == _LESS_THAN ? _LESS_AND_GREATER_THAN : _GREATER_THAN
        info.lower_bound_if_bounded = s.lower
    elseif S <: MOI.EqualTo{Float64}
        _throw_if_existing_lower(info.bound, info.type, S, f)
        _throw_if_existing_upper(info.bound, info.type, S, f)
        info.bound = _EQUAL_TO
        info.upper_bound_if_bounded = info.lower_bound_if_bounded = s.value
    else
        @assert S <: MOI.Interval{Float64}
        _throw_if_existing_lower(info.bound, info.type, S, f)
        _throw_if_existing_upper(info.bound, info.type, S, f)
        info.bound = _INTERVAL
        info.upper_bound_if_bounded = s.upper
        info.lower_bound_if_bounded = s.lower
    end
    index = MOI.ConstraintIndex{MOI.VariableIndex,typeof(s)}(f.value)
    # This sets the bounds in the inner model and set the cache in _VariableInfo
    # again (we could just set them there, but then _VariableInfo is in a
    # invalid state that trigger some asserts, i.e., has bound but no cache).
    MOI.set(model, MOI.ConstraintSet(), index, s)
    return index
end

function MOI.add_constraints(
    model::Optimizer,
    f::Vector{MOI.VariableIndex},
    s::Vector{S},
) where {S<:_SCALAR_SETS}
    for (fi, si) in zip(f, s)
        info = _info(model, fi)
        if S <: MOI.LessThan{Float64}
            _throw_if_existing_upper(info.bound, info.type, S, fi)
            info.bound =
                info.bound == _GREATER_THAN ? _LESS_AND_GREATER_THAN :
                _LESS_THAN
            info.upper_bound_if_bounded = si.upper
        elseif S <: MOI.GreaterThan{Float64}
            _throw_if_existing_lower(info.bound, info.type, S, fi)
            info.bound =
                info.bound == _LESS_THAN ? _LESS_AND_GREATER_THAN :
                _GREATER_THAN
            info.lower_bound_if_bounded = si.lower
        elseif S <: MOI.EqualTo{Float64}
            _throw_if_existing_lower(info.bound, info.type, S, fi)
            _throw_if_existing_upper(info.bound, info.type, S, fi)
            info.bound = _EQUAL_TO
            info.upper_bound_if_bounded = info.lower_bound_if_bounded = si.value
        else
            @assert S <: MOI.Interval{Float64}
            _throw_if_existing_lower(info.bound, info.type, S, fi)
            _throw_if_existing_upper(info.bound, info.type, S, fi)
            info.bound = _INTERVAL
            info.upper_bound_if_bounded = si.upper
            info.lower_bound_if_bounded = si.lower
        end
    end
    indices =
        [MOI.ConstraintIndex{MOI.VariableIndex,eltype(s)}(fi.value) for fi in f]
    _set_bounds(model, indices, s)
    return indices
end

function MOI.delete(
    model::Optimizer,
    c::MOI.ConstraintIndex{MOI.VariableIndex,MOI.LessThan{Float64}},
)
    MOI.throw_if_not_valid(model, c)
    info = _info(model, c)
    ret = GRBsetdblattrelement(model, "UB", Cint(info.column - 1), Inf)
    _check_ret(model, ret)
    _require_update(model)
    if info.bound == _LESS_AND_GREATER_THAN
        info.bound = _GREATER_THAN
    else
        info.bound = _NONE
    end
    info.upper_bound_if_bounded = NaN
    info.lessthan_name = ""
    model.name_to_constraint_index = nothing
    return
end

"""
    _set_variable_lower_bound(model, info, value)

This function is used to indirectly set the lower bound of a variable.

We need to do it this way to account for potential lower bounds of 0.0 added by
VectorOfVariables-in-SecondOrderCone constraints.

This does not look at `info.bound` and does not update
`info.lower_bound_if_bounded`.

See also `_get_variable_lower_bound`.
"""
function _set_variable_lower_bound(model, info, value)
    if info.num_soc_constraints == 0
        # No SOC constraints, set directly.
        @assert isnan(info.lower_bound_if_soc)
        GRBsetdblattrelement(model, "LB", Cint(info.column - 1), value)
        _require_update(model)
    elseif value >= 0.0
        # Regardless of whether there are SOC constraints, this is a valid
        # bound for the SOC constraint and should override any previous bounds.
        info.lower_bound_if_soc = NaN
        GRBsetdblattrelement(model, "LB", Cint(info.column - 1), value)
        _require_update(model)
    elseif isnan(info.lower_bound_if_soc)
        # Previously, we had a non-negative lower bound (i.e., it was set in
        # the case above). Now we are setting this with a negative one, but
        # there are still some SOC constraints, so we cache `value` and set the
        # variable lower bound to `0.0`.
        @assert value < 0.0
        GRBsetdblattrelement(model, "LB", Cint(info.column - 1), 0.0)
        _require_update(model)
        info.lower_bound_if_soc = value
    else
        # Previously, we had a negative lower bound. We are setting this with
        # another negative one, but there are still some SOC constraints.
        @assert info.lower_bound_if_soc < 0.0
        info.lower_bound_if_soc = value
    end
    return
end

"""
    _get_variable_lower_bound(model, info)

Get the current variable lower bound, ignoring a potential bound of `0.0` set
by a second order cone constraint, if an adequate `VariableIndex` constraint
is set (i.e., `info.bound` is not `_NONE` or `_LESS_THAN`) then use a cached
value; otherwise update the model if necessary and query the LB from it.

See also `_set_variable_lower_bound`.
"""
function _get_variable_lower_bound(model, info)
    if !isnan(info.lower_bound_if_soc)
        # There is a value stored. That means that we must have set a value
        # that was < 0.
        @assert info.lower_bound_if_soc < 0.0
        return info.lower_bound_if_soc
    elseif !isnan(info.lower_bound_if_bounded)
        @assert info.bound in
                (_GREATER_THAN, _LESS_AND_GREATER_THAN, _EQUAL_TO, _INTERVAL)
        return info.lower_bound_if_bounded
    end
    _update_if_necessary(model)
    valueP = Ref{Cdouble}()
    ret = GRBgetdblattrelement(model, "LB", Cint(info.column - 1), valueP)
    _check_ret(model, ret)
    return valueP[]
end

"""
    _get_variable_upper_bound(model, info)

Get the current variable upper bound, if an adequate `VariableIndex`
constraint is set (i.e., `info.bound` is not `_NONE` or `_GREATER_THAN`) then use
a cached value; otherwise update the model if necessary and query the UB from
it.

See also `_get_variable_lower_bound`.
"""
function _get_variable_upper_bound(model, info)
    if !isnan(info.upper_bound_if_bounded)
        @assert info.bound in
                (_LESS_THAN, _LESS_AND_GREATER_THAN, _EQUAL_TO, _INTERVAL)
        return info.upper_bound_if_bounded
    end
    _update_if_necessary(model)
    valueP = Ref{Cdouble}()
    ret = GRBgetdblattrelement(model, "UB", Cint(info.column - 1), valueP)
    _check_ret(model, ret)
    return valueP[]
end

function MOI.delete(
    model::Optimizer,
    c::MOI.ConstraintIndex{MOI.VariableIndex,MOI.GreaterThan{Float64}},
)
    MOI.throw_if_not_valid(model, c)
    info = _info(model, c)
    _set_variable_lower_bound(model, info, -Inf)
    if info.bound == _LESS_AND_GREATER_THAN
        info.bound = _LESS_THAN
    else
        info.bound = _NONE
    end
    info.lower_bound_if_bounded = NaN
    info.greaterthan_interval_or_equalto_name = ""
    model.name_to_constraint_index = nothing
    return
end

function MOI.delete(
    model::Optimizer,
    c::MOI.ConstraintIndex{MOI.VariableIndex,MOI.Interval{Float64}},
)
    MOI.throw_if_not_valid(model, c)
    info = _info(model, c)
    _set_variable_lower_bound(model, info, -Inf)
    ret = GRBsetdblattrelement(model, "UB", Cint(info.column - 1), Inf)
    _check_ret(model, ret)
    _require_update(model)
    info.bound = _NONE
    info.upper_bound_if_bounded = info.lower_bound_if_bounded = NaN
    info.greaterthan_interval_or_equalto_name = ""
    model.name_to_constraint_index = nothing
    return
end

function MOI.delete(
    model::Optimizer,
    c::MOI.ConstraintIndex{MOI.VariableIndex,MOI.EqualTo{Float64}},
)
    MOI.throw_if_not_valid(model, c)
    info = _info(model, c)
    _set_variable_lower_bound(model, info, -Inf)
    ret = GRBsetdblattrelement(model, "UB", Cint(info.column - 1), Inf)
    _check_ret(model, ret)
    _require_update(model)
    info.bound = _NONE
    info.upper_bound_if_bounded = info.lower_bound_if_bounded = NaN
    info.greaterthan_interval_or_equalto_name = ""
    model.name_to_constraint_index = nothing
    return
end

function MOI.get(
    model::Optimizer,
    ::MOI.ConstraintSet,
    c::MOI.ConstraintIndex{MOI.VariableIndex,MOI.GreaterThan{Float64}},
)
    MOI.throw_if_not_valid(model, c)
    info = _info(model, c)
    @assert !isnan(info.lower_bound_if_bounded)
    return MOI.GreaterThan(_get_variable_lower_bound(model, info))
end

function MOI.get(
    model::Optimizer,
    ::MOI.ConstraintSet,
    c::MOI.ConstraintIndex{MOI.VariableIndex,MOI.LessThan{Float64}},
)
    MOI.throw_if_not_valid(model, c)
    info = _info(model, c)
    @assert !isnan(info.upper_bound_if_bounded)
    return MOI.LessThan(_get_variable_upper_bound(model, info))
end

function MOI.get(
    model::Optimizer,
    ::MOI.ConstraintSet,
    c::MOI.ConstraintIndex{MOI.VariableIndex,MOI.EqualTo{Float64}},
)
    MOI.throw_if_not_valid(model, c)
    info = _info(model, c)
    @assert !isnan(info.upper_bound_if_bounded)
    @assert info.upper_bound_if_bounded == info.lower_bound_if_bounded
    return MOI.EqualTo(_get_variable_lower_bound(model, info))
end

function MOI.get(
    model::Optimizer,
    ::MOI.ConstraintSet,
    c::MOI.ConstraintIndex{MOI.VariableIndex,MOI.Interval{Float64}},
)
    MOI.throw_if_not_valid(model, c)
    info = _info(model, c)
    @assert !isnan(info.upper_bound_if_bounded)
    @assert !isnan(info.lower_bound_if_bounded)
    return MOI.Interval(
        _get_variable_lower_bound(model, info),
        _get_variable_upper_bound(model, info),
    )
end

function _set_bounds(
    model::Optimizer,
    indices::Vector{MOI.ConstraintIndex{MOI.VariableIndex,S}},
    sets::Vector{S},
) where {S}
    lower_columns, lower_values = Cint[], Float64[]
    upper_columns, upper_values = Cint[], Float64[]
    for (c, s) in zip(indices, sets)
        lower, upper = _bounds(s)
        info = _info(model, c)
        if lower !== nothing
            if info.num_soc_constraints == 0
                push!(lower_columns, Cint(info.column - 1))
                push!(lower_values, lower)
            else
                _set_variable_lower_bound(model, info, lower)
            end
        end
        if upper !== nothing
            push!(upper_columns, Cint(info.column - 1))
            push!(upper_values, upper)
        end
    end
    if length(lower_columns) > 0
        ret = GRBsetdblattrlist(
            model,
            "LB",
            length(lower_columns),
            lower_columns,
            lower_values,
        )
        _check_ret(model, ret)
    end
    if length(upper_columns) > 0
        ret = GRBsetdblattrlist(
            model,
            "UB",
            length(upper_columns),
            upper_columns,
            upper_values,
        )
        _check_ret(model, ret)
    end
    return _require_update(model)
end

function MOI.set(
    model::Optimizer,
    ::MOI.ConstraintSet,
    c::MOI.ConstraintIndex{MOI.VariableIndex,S},
    s::S,
) where {S<:_SCALAR_SETS}
    MOI.throw_if_not_valid(model, c)
    lower, upper = _bounds(s)
    info = _info(model, c)
    if lower !== nothing
        @assert !isnan(info.lower_bound_if_bounded)
        info.lower_bound_if_bounded = lower
        _set_variable_lower_bound(model, info, lower)
    end
    if upper !== nothing
        @assert !isnan(info.upper_bound_if_bounded)
        info.upper_bound_if_bounded = upper
        ret = GRBsetdblattrelement(model, "UB", Cint(info.column - 1), upper)
        _check_ret(model, ret)
    end
    _require_update(model)
    return
end

function MOI.add_constraint(
    model::Optimizer,
    f::MOI.VariableIndex,
    ::MOI.ZeroOne,
)
    info = _info(model, f)
    ret =
        GRBsetcharattrelement(model, "VType", Cint(info.column - 1), Char('B'))
    _check_ret(model, ret)
    _require_update(model)
    info.type = GRB_BINARY
    return MOI.ConstraintIndex{MOI.VariableIndex,MOI.ZeroOne}(f.value)
end

function MOI.delete(
    model::Optimizer,
    c::MOI.ConstraintIndex{MOI.VariableIndex,MOI.ZeroOne},
)
    MOI.throw_if_not_valid(model, c)
    info = _info(model, c)
    ret =
        GRBsetcharattrelement(model, "VType", Cint(info.column - 1), Char('C'))
    _check_ret(model, ret)
    _require_update(model)
    info.type = GRB_CONTINUOUS
    info.type_constraint_name = ""
    model.name_to_constraint_index = nothing
    return
end

function MOI.get(
    model::Optimizer,
    ::MOI.ConstraintSet,
    c::MOI.ConstraintIndex{MOI.VariableIndex,MOI.ZeroOne},
)
    MOI.throw_if_not_valid(model, c)
    return MOI.ZeroOne()
end

function MOI.add_constraint(
    model::Optimizer,
    f::MOI.VariableIndex,
    ::MOI.Integer,
)
    info = _info(model, f)
    ret =
        GRBsetcharattrelement(model, "VType", Cint(info.column - 1), Char('I'))
    _check_ret(model, ret)
    _require_update(model)
    info.type = GRB_INTEGER
    return MOI.ConstraintIndex{MOI.VariableIndex,MOI.Integer}(f.value)
end

function MOI.delete(
    model::Optimizer,
    c::MOI.ConstraintIndex{MOI.VariableIndex,MOI.Integer},
)
    MOI.throw_if_not_valid(model, c)
    info = _info(model, c)
    ret =
        GRBsetcharattrelement(model, "VType", Cint(info.column - 1), Char('C'))
    _check_ret(model, ret)
    _require_update(model)
    info.type = GRB_CONTINUOUS
    info.type_constraint_name = ""
    model.name_to_constraint_index = nothing
    return
end

function MOI.get(
    model::Optimizer,
    ::MOI.ConstraintSet,
    c::MOI.ConstraintIndex{MOI.VariableIndex,MOI.Integer},
)
    MOI.throw_if_not_valid(model, c)
    return MOI.Integer()
end

function MOI.add_constraint(
    model::Optimizer,
    f::MOI.VariableIndex,
    s::MOI.Semicontinuous{Float64},
)
    info = _info(model, f)
    _throw_if_existing_lower(info.bound, info.type, typeof(s), f)
    _throw_if_existing_upper(info.bound, info.type, typeof(s), f)
    ret =
        GRBsetcharattrelement(model, "VType", Cint(info.column - 1), Char('S'))
    _check_ret(model, ret)
    _set_variable_lower_bound(model, info, s.lower)
    ret = GRBsetdblattrelement(model, "UB", Cint(info.column - 1), s.upper)
    _check_ret(model, ret)
    _require_update(model)
    info.type = GRB_SEMICONT
    return MOI.ConstraintIndex{MOI.VariableIndex,MOI.Semicontinuous{Float64}}(
        f.value,
    )
end

function MOI.delete(
    model::Optimizer,
    c::MOI.ConstraintIndex{MOI.VariableIndex,MOI.Semicontinuous{Float64}},
)
    MOI.throw_if_not_valid(model, c)
    info = _info(model, c)
    ret =
        GRBsetcharattrelement(model, "VType", Cint(info.column - 1), Char('C'))
    _check_ret(model, ret)
    _set_variable_lower_bound(model, info, -Inf)
    ret = GRBsetdblattrelement(model, "UB", Cint(info.column - 1), Inf)
    _check_ret(model, ret)
    _require_update(model)
    info.type = GRB_CONTINUOUS
    info.type_constraint_name = ""
    model.name_to_constraint_index = nothing
    return
end

function MOI.get(
    model::Optimizer,
    ::MOI.ConstraintSet,
    c::MOI.ConstraintIndex{MOI.VariableIndex,MOI.Semicontinuous{Float64}},
)
    MOI.throw_if_not_valid(model, c)
    info = _info(model, c)
    _update_if_necessary(model)
    lower = _get_variable_lower_bound(model, info)
    upper = Ref{Cdouble}()
    ret = GRBgetdblattrelement(model, "UB", Cint(info.column - 1), upper)
    return MOI.Semicontinuous(lower, upper[])
end

function MOI.add_constraint(
    model::Optimizer,
    f::MOI.VariableIndex,
    s::MOI.Semiinteger{Float64},
)
    info = _info(model, f)
    _throw_if_existing_lower(info.bound, info.type, typeof(s), f)
    _throw_if_existing_upper(info.bound, info.type, typeof(s), f)
    ret =
        GRBsetcharattrelement(model, "VType", Cint(info.column - 1), Char('N'))
    _check_ret(model, ret)
    _set_variable_lower_bound(model, info, s.lower)
    ret = GRBsetdblattrelement(model, "UB", Cint(info.column - 1), s.upper)
    _check_ret(model, ret)
    _require_update(model)
    info.type = GRB_SEMIINT
    return MOI.ConstraintIndex{MOI.VariableIndex,MOI.Semiinteger{Float64}}(
        f.value,
    )
end

function MOI.delete(
    model::Optimizer,
    c::MOI.ConstraintIndex{MOI.VariableIndex,MOI.Semiinteger{Float64}},
)
    MOI.throw_if_not_valid(model, c)
    info = _info(model, c)
    ret =
        GRBsetcharattrelement(model, "VType", Cint(info.column - 1), Char('C'))
    _check_ret(model, ret)
    _set_variable_lower_bound(model, info, -Inf)
    ret = GRBsetdblattrelement(model, "UB", Cint(info.column - 1), Inf)
    _check_ret(model, ret)
    _require_update(model)
    info.type = GRB_CONTINUOUS
    info.type_constraint_name = ""
    model.name_to_constraint_index = nothing
    return
end

function MOI.get(
    model::Optimizer,
    ::MOI.ConstraintSet,
    c::MOI.ConstraintIndex{MOI.VariableIndex,MOI.Semiinteger{Float64}},
)
    MOI.throw_if_not_valid(model, c)
    info = _info(model, c)
    _update_if_necessary(model)
    lower = _get_variable_lower_bound(model, info)
    upper = Ref{Cdouble}()
    ret = GRBgetdblattrelement(model, "UB", Cint(info.column - 1), upper)
    return MOI.Semiinteger(lower, upper[])
end

function MOI.get(
    model::Optimizer,
    ::MOI.ConstraintName,
    c::MOI.ConstraintIndex{MOI.VariableIndex,S},
) where {S}
    MOI.throw_if_not_valid(model, c)
    info = _info(model, c)
    if S <: MOI.LessThan
        return info.lessthan_name
    elseif S <: Union{MOI.GreaterThan,MOI.Interval,MOI.EqualTo}
        return info.greaterthan_interval_or_equalto_name
    else
        @assert S <: Union{
            MOI.ZeroOne,
            MOI.Integer,
            MOI.Semiinteger,
            MOI.Semicontinuous,
        }
        return info.type_constraint_name
    end
end

function MOI.set(
    model::Optimizer,
    ::MOI.ConstraintName,
    c::MOI.ConstraintIndex{MOI.VariableIndex,S},
    name::String,
) where {S}
    MOI.throw_if_not_valid(model, c)
    info = _info(model, c)
    old_name = ""
    if S <: MOI.LessThan
        old_name = info.lessthan_name
        info.lessthan_name = name
    elseif S <: Union{MOI.GreaterThan,MOI.Interval,MOI.EqualTo}
        old_name = info.greaterthan_interval_or_equalto_name
        info.greaterthan_interval_or_equalto_name = name
    else
        @assert S <: Union{
            MOI.ZeroOne,
            MOI.Integer,
            MOI.Semiinteger,
            MOI.Semicontinuous,
        }
        info.type_constraint_name
        info.type_constraint_name = name
    end
    model.name_to_constraint_index = nothing
    return
end

###
### ScalarAffineFunction-in-Set
###

function _info(
    model::Optimizer,
    key::MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64},<:Any},
)
    if haskey(model.affine_constraint_info, key.value)
        return model.affine_constraint_info[key.value]
    end
    return throw(MOI.InvalidIndex(key))
end

function MOI.is_valid(
    model::Optimizer,
    c::MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64},S},
) where {S}
    info = get(model.affine_constraint_info, c.value, nothing)
    if info === nothing
        return false
    else
        return typeof(info.set) == S
    end
end

function MOI.add_constraint(
    model::Optimizer,
    f::MOI.ScalarAffineFunction{Float64},
    s::_SUPPORTED_SCALAR_SETS,
)
    if !iszero(f.constant)
        throw(
            MOI.ScalarFunctionConstantNotZero{Float64,typeof(f),typeof(s)}(
                f.constant,
            ),
        )
    end
    model.last_constraint_index += 1
    model.affine_constraint_info[model.last_constraint_index] =
        _ConstraintInfo(length(model.affine_constraint_info) + 1, s)

    indices, coefficients = _indices_and_coefficients(model, f)
    sense, rhs = _sense_and_rhs(s)
    ret = GRBaddconstr(
        model,
        length(indices),
        indices,
        coefficients,
        sense,
        rhs,
        "",
    )
    _check_ret(model, ret)
    _require_update(model)
    return MOI.ConstraintIndex{typeof(f),typeof(s)}(model.last_constraint_index)
end

function MOI.add_constraints(
    model::Optimizer,
    f::Vector{MOI.ScalarAffineFunction{Float64}},
    s::Vector{<:_SUPPORTED_SCALAR_SETS},
)
    if length(f) != length(s)
        error("Number of functions does not equal number of sets.")
    end
    canonicalized_functions = MOI.Utilities.canonical.(f)
    # First pass: compute number of non-zeros to allocate space.
    nnz = 0
    for fi in canonicalized_functions
        if !iszero(fi.constant)
            throw(
                MOI.ScalarFunctionConstantNotZero{Float64,eltype(f),eltype(s)}(
                    fi.constant,
                ),
            )
        end
        nnz += length(fi.terms)
    end
    # Initialize storage
    indices = Vector{MOI.ConstraintIndex{eltype(f),eltype(s)}}(undef, length(f))
    row_starts = Vector{Cint}(undef, length(f) + 1)
    row_starts[1] = 0
    columns = Vector{Cint}(undef, nnz)
    coefficients = Vector{Float64}(undef, nnz)
    senses = Vector{Cchar}(undef, length(f))
    rhss = Vector{Float64}(undef, length(f))
    # Second pass: loop through, passing views to _indices_and_coefficients.
    for (i, (fi, si)) in enumerate(zip(canonicalized_functions, s))
        senses[i], rhss[i] = _sense_and_rhs(si)
        row_starts[i+1] = row_starts[i] + length(fi.terms)
        _indices_and_coefficients(
            view(columns, (1+row_starts[i]):row_starts[i+1]),
            view(coefficients, (1+row_starts[i]):row_starts[i+1]),
            model,
            fi,
        )
        model.last_constraint_index += 1
        indices[i] = MOI.ConstraintIndex{eltype(f),eltype(s)}(
            model.last_constraint_index,
        )
        model.affine_constraint_info[model.last_constraint_index] =
            _ConstraintInfo(length(model.affine_constraint_info) + 1, si)
    end
    pop!(row_starts)  # Gurobi doesn't need the final row start.
    ret = GC.@preserve GRBaddconstrs(
        model,
        length(f),
        length(coefficients),
        row_starts,
        columns,
        coefficients,
        senses,
        rhss,
        C_NULL,
    )
    _require_update(model)
    return indices
end

function MOI.delete(
    model::Optimizer,
    cs::Vector{<:MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64},<:Any}},
)
    # The `sort!`s called are necessary for improving the efficiency of the
    # updates (i.e., without them we cannot use `searchsorted*` methods).
    _update_if_necessary(model)
    rows_to_delete = sort!([Cint(_info(model, x).row - 1) for x in cs])
    ret = GRBdelconstrs(model, length(rows_to_delete), rows_to_delete)
    _require_update(model)
    for (_, info) in model.affine_constraint_info
        # The trick here is: searchsortedlast returns, in O(log n), the
        # last index with a row smaller than info.row, over rows_to_delete
        # this is the same as the number of rows deleted before it, and
        # how much its value need to be shifted.
        info.row -= searchsortedlast(rows_to_delete, info.row - 1)
    end
    cs_values = sort!(getfield.(cs, :value))
    # If the key of an model.affine_constraint_info entry is in cs_values,
    # then that entry is deleted.
    filter!(model.affine_constraint_info) do pair
        return isempty(searchsorted(cs_values, pair.first))
    end
    model.name_to_constraint_index = nothing
    return
end

function MOI.delete(
    model::Optimizer,
    c::MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64},<:Any},
)
    row = _info(model, c).row
    _update_if_necessary(model)
    ret = GRBdelconstrs(model, 1, Ref{Cint}(row - 1))
    _check_ret(model, ret)
    _require_update(model)
    for (key, info) in model.affine_constraint_info
        if info.row > row
            info.row -= 1
        end
    end
    delete!(model.affine_constraint_info, c.value)
    model.name_to_constraint_index = nothing
    return
end

function MOI.get(
    model::Optimizer,
    ::MOI.ConstraintSet,
    c::MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64},S},
) where {S}
    _update_if_necessary(model)
    rhs = Ref{Cdouble}()
    ret = GRBgetdblattrelement(model, "RHS", Cint(_info(model, c).row - 1), rhs)
    _check_ret(model, ret)
    return S(rhs[])
end

function MOI.set(
    model::Optimizer,
    ::MOI.ConstraintSet,
    c::MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64},S},
    s::S,
) where {S}
    ret = GRBsetdblattrelement(
        model,
        "RHS",
        Cint(_info(model, c).row - 1),
        MOI.constant(s),
    )
    _check_ret(model, ret)
    _require_update(model)
    return
end

function MOI.get(
    model::Optimizer,
    ::MOI.ConstraintFunction,
    c::MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64},S},
) where {S}
    _update_if_necessary(model)
    numnzP = Ref{Cint}()
    row = Cint(_info(model, c).row - 1)
    ret = GRBgetconstrs(model, numnzP, C_NULL, C_NULL, C_NULL, row, 1)
    _check_ret(model, ret)
    cbeg = Array{Cint}(undef, 2)
    cind = Array{Cint}(undef, numnzP[])
    cval = Array{Cdouble}(undef, numnzP[])
    ret = GRBgetconstrs(model, numnzP, cbeg, cind, cval, row, 1)
    _check_ret(model, ret)
    terms = Vector{MOI.ScalarAffineTerm{Float64}}(undef, length(cind))
    for i in 1:length(terms)
        terms[i] = MOI.ScalarAffineTerm(
            cval[i],
            model.variable_info[CleverDicts.LinearIndex(cind[i] + 1)].index,
        )
    end
    return MOI.ScalarAffineFunction(terms, 0.0)
end

function MOI.get(
    model::Optimizer,
    ::MOI.ConstraintName,
    c::MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64},<:Any},
)
    return _info(model, c).name
end

function MOI.set(
    model::Optimizer,
    ::MOI.ConstraintName,
    c::MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64},<:Any},
    name::String,
)
    info = _info(model, c)
    info.name = name
    if !isempty(name)
        ret =
            GRBsetstrattrelement(model, "ConstrName", Cint(info.row - 1), name)
        _check_ret(model, ret)
        _require_update(model)
    end
    model.name_to_constraint_index = nothing
    return
end

function MOI.get(model::Optimizer, ::Type{MOI.ConstraintIndex}, name::String)
    if model.name_to_constraint_index === nothing
        _rebuild_name_to_constraint_index(model)
    end
    if haskey(model.name_to_constraint_index, name)
        constr = model.name_to_constraint_index[name]
        if constr === nothing
            error("Duplicate constraint name detected: $(name)")
        end
        return constr
    end
    return nothing
end

function MOI.get(
    model::Optimizer,
    C::Type{MOI.ConstraintIndex{F,S}},
    name::String,
) where {F,S}
    index = MOI.get(model, MOI.ConstraintIndex, name)
    if typeof(index) == C
        return index::MOI.ConstraintIndex{F,S}
    end
    return nothing
end

function _rebuild_name_to_constraint_index(model::Optimizer)
    model.name_to_constraint_index =
        Dict{String,Union{Nothing,MOI.ConstraintIndex}}()
    _rebuild_name_to_constraint_index_util(
        model,
        model.affine_constraint_info,
        MOI.ScalarAffineFunction{Float64},
    )
    _rebuild_name_to_constraint_index_util(
        model,
        model.quadratic_constraint_info,
        MOI.ScalarQuadraticFunction{Float64},
    )
    _rebuild_name_to_constraint_index_util(
        model,
        model.sos_constraint_info,
        MOI.VectorOfVariables,
    )
    _rebuild_name_to_constraint_index_variables(model)
    return
end

function _rebuild_name_to_constraint_index_util(model::Optimizer, dict, F)
    for (index, info) in dict
        if info.name == ""
            continue
        elseif haskey(model.name_to_constraint_index, info.name)
            model.name_to_constraint_index[info.name] = nothing
        else
            model.name_to_constraint_index[info.name] =
                MOI.ConstraintIndex{F,typeof(info.set)}(index)
        end
    end
    return
end

function _rebuild_name_to_constraint_index_variables(model::Optimizer)
    for (key, info) in model.variable_info
        for S in (
            MOI.LessThan{Float64},
            MOI.GreaterThan{Float64},
            MOI.EqualTo{Float64},
            MOI.Interval{Float64},
            MOI.ZeroOne,
            MOI.Integer,
            MOI.Semicontinuous{Float64},
            MOI.Semiinteger{Float64},
        )
            constraint_name = ""
            if info.bound in _bound_enums(S)
                constraint_name =
                    S == MOI.LessThan{Float64} ? info.lessthan_name :
                    info.greaterthan_interval_or_equalto_name
            elseif info.type in _type_enums(S)
                constraint_name = info.type_constraint_name
            end
            if constraint_name == ""
                continue
            elseif haskey(model.name_to_constraint_index, constraint_name)
                model.name_to_constraint_index[constraint_name] = nothing
            else
                model.name_to_constraint_index[constraint_name] =
                    MOI.ConstraintIndex{MOI.VariableIndex,S}(key.value)
            end
        end
    end
    return
end

###
### ScalarQuadraticFunction-in-SCALAR_SET
###

function _info(
    model::Optimizer,
    c::MOI.ConstraintIndex{MOI.ScalarQuadraticFunction{Float64},S},
) where {S}
    if haskey(model.quadratic_constraint_info, c.value)
        return model.quadratic_constraint_info[c.value]
    end
    return throw(MOI.InvalidIndex(c))
end

function MOI.add_constraint(
    model::Optimizer,
    f::MOI.ScalarQuadraticFunction{Float64},
    s::_SCALAR_SETS,
)
    if !iszero(f.constant)
        throw(
            MOI.ScalarFunctionConstantNotZero{Float64,typeof(f),typeof(s)}(
                f.constant,
            ),
        )
    end
    indices, coefficients, I, J, V = _indices_and_coefficients(model, f)
    sense, rhs = _sense_and_rhs(s)
    ret = GRBaddqconstr(
        model,
        length(indices),
        indices,
        coefficients,
        length(I),
        I,
        J,
        V,
        sense,
        rhs,
        "",
    )
    _check_ret(model, ret)
    _require_update(model)
    model.last_constraint_index += 1
    model.quadratic_constraint_info[model.last_constraint_index] =
        _ConstraintInfo(length(model.quadratic_constraint_info) + 1, s)
    return MOI.ConstraintIndex{MOI.ScalarQuadraticFunction{Float64},typeof(s)}(
        model.last_constraint_index,
    )
end

function MOI.is_valid(
    model::Optimizer,
    c::MOI.ConstraintIndex{MOI.ScalarQuadraticFunction{Float64},S},
) where {S}
    info = get(model.quadratic_constraint_info, c.value, nothing)
    return info !== nothing && typeof(info.set) == S
end

function MOI.delete(
    model::Optimizer,
    c::MOI.ConstraintIndex{MOI.ScalarQuadraticFunction{Float64},S},
) where {S}
    _update_if_necessary(model)
    info = _info(model, c)
    ret = GRBdelqconstrs(model, 1, Ref{Cint}(info.row - 1))
    _check_ret(model, ret)
    _require_update(model)
    for (key, info_2) in model.quadratic_constraint_info
        if info_2.row > info.row
            info_2.row -= 1
        end
    end
    delete!(model.quadratic_constraint_info, c.value)
    model.name_to_constraint_index = nothing
    return
end

function MOI.get(
    model::Optimizer,
    ::MOI.ConstraintSet,
    c::MOI.ConstraintIndex{MOI.ScalarQuadraticFunction{Float64},S},
) where {S}
    _update_if_necessary(model)
    rhs = Ref{Cdouble}()
    ret =
        GRBgetdblattrelement(model, "QCRHS", Cint(_info(model, c).row - 1), rhs)
    _check_ret(model, ret)
    return S(rhs[])
end

function _getqconstr(model::Optimizer, row::Cint)
    numlnzP = Ref{Cint}()
    numqnzP = Ref{Cint}()
    ret = GRBgetqconstr(
        model,
        row,
        numlnzP,
        C_NULL,
        C_NULL,
        numqnzP,
        C_NULL,
        C_NULL,
        C_NULL,
    )
    _check_ret(model, ret)
    lind = Array{Cint}(undef, numlnzP[])
    lval = Array{Cdouble}(undef, numlnzP[])
    qrow = Array{Cint}(undef, numqnzP[])
    qcol = Array{Cint}(undef, numqnzP[])
    qval = Array{Cdouble}(undef, numqnzP[])
    ret = GRBgetqconstr(
        model,
        row,
        numlnzP,
        lind,
        lval,
        numqnzP,
        qrow,
        qcol,
        qval,
    )
    _check_ret(model, ret)
    return lind, lval, qrow, qcol, qval
end

function MOI.get(
    model::Optimizer,
    ::MOI.ConstraintFunction,
    c::MOI.ConstraintIndex{MOI.ScalarQuadraticFunction{Float64},S},
) where {S}
    _update_if_necessary(model)
    li, lv, I, J, V = _getqconstr(model, Cint(_info(model, c).row - 1))
    affine_terms = Vector{MOI.ScalarAffineTerm{Float64}}(undef, length(li))
    for i in 1:length(li)
        affine_terms[i] = MOI.ScalarAffineTerm(
            lv[i],
            model.variable_info[CleverDicts.LinearIndex(li[i] + 1)].index,
        )
    end
    quadratic_terms = Vector{MOI.ScalarQuadraticTerm{Float64}}(undef, length(I))
    for i in 1:length(I)
        v = I[i] == J[i] ? 2V[i] : V[i]
        quadratic_terms[i] = MOI.ScalarQuadraticTerm(
            v,
            model.variable_info[CleverDicts.LinearIndex(I[i] + 1)].index,
            model.variable_info[CleverDicts.LinearIndex(J[i] + 1)].index,
        )
    end
    return MOI.ScalarQuadraticFunction(quadratic_terms, affine_terms, 0.0)
end

function MOI.get(
    model::Optimizer,
    ::MOI.ConstraintName,
    c::MOI.ConstraintIndex{MOI.ScalarQuadraticFunction{Float64},S},
) where {S}
    return _info(model, c).name
end

function MOI.set(
    model::Optimizer,
    ::MOI.ConstraintName,
    c::MOI.ConstraintIndex{MOI.ScalarQuadraticFunction{Float64},S},
    name::String,
) where {S}
    info = _info(model, c)
    info.name = name
    _update_if_necessary(model)
    ret = GRBsetstrattrelement(model, "QCName", Cint(info.row - 1), name)
    _check_ret(model, ret)
    _require_update(model)
    model.name_to_constraint_index = nothing
    return
end

###
### VectorOfVariables-in-SOS{I|II}
###

const _SOS = Union{MOI.SOS1{Float64},MOI.SOS2{Float64}}

function _info(
    model::Optimizer,
    key::MOI.ConstraintIndex{MOI.VectorOfVariables,<:_SOS},
)
    if haskey(model.sos_constraint_info, key.value)
        return model.sos_constraint_info[key.value]
    end
    return throw(MOI.InvalidIndex(key))
end

_sos_type(::MOI.SOS1) = GRB_SOS_TYPE1
_sos_type(::MOI.SOS2) = GRB_SOS_TYPE2

function MOI.is_valid(
    model::Optimizer,
    c::MOI.ConstraintIndex{MOI.VectorOfVariables,S},
) where {S}
    info = get(model.sos_constraint_info, c.value, nothing)
    if info === nothing || typeof(info.set) != S
        return false
    end
    f = MOI.get(model, MOI.ConstraintFunction(), c)
    return all(MOI.is_valid.(model, f.variables))
end

function MOI.add_constraint(model::Optimizer, f::MOI.VectorOfVariables, s::_SOS)
    columns = Cint[column(model, v) - 1 for v in f.variables]
    ret = GRBaddsos(
        model,
        1,
        length(columns),
        Ref{Cint}(_sos_type(s)),
        Ref{Cint}(0),
        columns,
        s.weights,
    )
    _check_ret(model, ret)
    model.last_constraint_index += 1
    index = MOI.ConstraintIndex{MOI.VectorOfVariables,typeof(s)}(
        model.last_constraint_index,
    )
    model.sos_constraint_info[index.value] =
        _ConstraintInfo(length(model.sos_constraint_info) + 1, s)
    _require_update(model)
    return index
end

function MOI.delete(
    model::Optimizer,
    c::MOI.ConstraintIndex{MOI.VectorOfVariables,<:_SOS},
)
    row = _info(model, c).row
    _update_if_necessary(model)
    rows = Ref{Cint}(row - 1)
    ret = GRBdelsos(model, 1, rows)
    _check_ret(model, ret)
    _require_update(model)
    for (key, info) in model.sos_constraint_info
        if info.row > row
            info.row -= 1
        end
    end
    delete!(model.sos_constraint_info, c.value)
    model.name_to_constraint_index = nothing
    return
end

function MOI.delete(
    model::Optimizer,
    cs::Vector{<:MOI.ConstraintIndex{MOI.VectorOfVariables,<:_SOS}},
)
    # The batch delete method is same as batch delete of `MOI.ScalarAffineFunction{Float64}, <:Any}`
    _update_if_necessary(model)
    rows_to_delete = sort!([Cint(_info(model, c).row - 1) for c in cs])
    ret = GRBdelsos(model, length(rows_to_delete), rows_to_delete)
    _check_ret(model, ret)
    _require_update(model)
    for (_, info) in model.sos_constraint_info
        info.row -= searchsortedlast(rows_to_delete, info.row - 1)
    end
    cs_values = sort!(getfield.(cs, :value))
    filter!(model.sos_constraint_info) do pair
        return isempty(searchsorted(cs_values, pair.first))
    end
    model.name_to_constraint_index = nothing
    return
end

function MOI.get(
    model::Optimizer,
    ::MOI.ConstraintName,
    c::MOI.ConstraintIndex{MOI.VectorOfVariables,<:Any},
)
    return _info(model, c).name
end

function MOI.set(
    model::Optimizer,
    ::MOI.ConstraintName,
    c::MOI.ConstraintIndex{MOI.VectorOfVariables,<:Any},
    name::String,
)
    info = _info(model, c)
    info.name = name
    model.name_to_constraint_index = nothing
    return
end

function MOI.get(
    model::Optimizer,
    ::MOI.ConstraintSet,
    c::MOI.ConstraintIndex{MOI.VectorOfVariables,S},
) where {S<:_SOS}
    _update_if_necessary(model)
    nummembersP = Ref{Cint}()
    start = Cint(_info(model, c).row - 1)
    sostype = Ref{Cint}()
    rbeg = Vector{Cint}(undef, 2)
    ret = GRBgetsos(model, nummembersP, sostype, rbeg, C_NULL, C_NULL, start, 1)
    _check_ret(model, ret)
    ind = Array{Cint}(undef, nummembersP[])
    weight = Array{Cdouble}(undef, nummembersP[])
    ret = GRBgetsos(model, nummembersP, sostype, rbeg, ind, weight, start, 1)
    _check_ret(model, ret)
    return S(weight)
end

function MOI.get(
    model::Optimizer,
    ::MOI.ConstraintFunction,
    c::MOI.ConstraintIndex{MOI.VectorOfVariables,S},
) where {S<:_SOS}
    _update_if_necessary(model)
    nummembersP = Ref{Cint}()
    start = Cint(_info(model, c).row - 1)
    sostype = Ref{Cint}()
    rbeg = Vector{Cint}(undef, 2)
    ret = GRBgetsos(model, nummembersP, sostype, rbeg, C_NULL, C_NULL, start, 1)
    _check_ret(model, ret)
    ind = Array{Cint}(undef, nummembersP[])
    weight = Array{Cdouble}(undef, nummembersP[])
    ret = GRBgetsos(model, nummembersP, sostype, rbeg, ind, weight, start, 1)
    _check_ret(model, ret)
    return MOI.VectorOfVariables([
        model.variable_info[CleverDicts.LinearIndex(i + 1)].index for i in ind
    ])
end

###
### Optimize methods.
###

function _check_moi_callback_validity(model::Optimizer)
    has_moi_callback =
        model.lazy_callback !== nothing ||
        model.user_cut_callback !== nothing ||
        model.heuristic_callback !== nothing
    if has_moi_callback && model.has_generic_callback
        error(
            "Cannot use Gurobi.CallbackFunction as well as " *
            "MOI.AbstractCallbackFunction",
        )
    end
    return has_moi_callback
end

function MOI.optimize!(model::Optimizer)
    # Note: although Gurobi will call update regardless, we do it now so that
    # the appropriate `needs_update` flag is set.
    _update_if_necessary(model)

    # Initialize callbacks if necessary.
    if _check_moi_callback_validity(model)
        MOI.set(model, CallbackFunction(), _default_moi_callback(model))
        model.has_generic_callback = false
    elseif model.enable_interrupts && !model.has_generic_callback
        # From the docstring of disable_sigint, "External functions that do not
        # call julia code or julia runtime automatically disable sigint during
        # their execution." We don't want this though! We want to be able to
        # SIGINT Gurobi, and then catch it as an interrupt. As a hack, until
        # Julia introduces an interruptible ccall --- which it likely won't
        # https://github.com/JuliaLang/julia/issues/2622 --- set a null
        # callback.
        MOI.set(model, CallbackFunction(), (x, y) -> nothing)
    end

    # Catch [CTRL+C], even when Julia is run from a script not in interactive
    # mode. If `true`, then a script would call `atexit` without throwing the
    # `InterruptException`. `false` is the default in interactive mode.
    #
    # TODO(odow): Julia 1.5 exposes `Base.exit_on_sigint(::Bool)`.
    ccall(:jl_exit_on_sigint, Cvoid, (Cint,), false)
    model.ret_GRBoptimize = GRBoptimize(model)
    _check_ret_GRBoptimize(model)
    if !isinteractive()
        ccall(:jl_exit_on_sigint, Cvoid, (Cint,), true)
    end

    # Post-optimize caching to speed up the checks in VariablePrimal and
    # ConstraintDual.
    model.has_infeasibility_cert =
        MOI.get(model, MOI.DualStatus()) == MOI.INFEASIBILITY_CERTIFICATE
    model.has_unbounded_ray =
        MOI.get(model, MOI.PrimalStatus()) == MOI.INFEASIBILITY_CERTIFICATE

    return
end

function _throw_if_optimize_in_progress(model, attr)
    if model.callback_state != _CB_NONE
        throw(MOI.OptimizeInProgress(attr))
    end
end

# These strings are taken directly from the following page of the online Gurobi
# documentation:
# https://www.com/documentation/8.1/refman/optimization_status_codes.html#sec:StatusCodes
const _RAW_STATUS_STRINGS = [
    # TerminationStatus,          RawStatusString
    (
        MOI.OPTIMIZE_NOT_CALLED,
        "Model is loaded, but no solution information is available.",
    ),
    (
        MOI.OPTIMAL,
        "Model was solved to optimality (subject to tolerances), and an optimal solution is available.",
    ),
    (MOI.INFEASIBLE, "Model was proven to be infeasible."),
    (
        MOI.INFEASIBLE_OR_UNBOUNDED,
        "Model was proven to be either infeasible or unbounded. To obtain a more definitive conclusion, set the DualReductions parameter to 0 and reoptimize.",
    ),
    (
        MOI.DUAL_INFEASIBLE,
        "Model was proven to be unbounded. Important note: an unbounded status indicates the presence of an unbounded ray that allows the objective to improve without limit. It says nothing about whether the model has a feasible solution. If you require information on feasibility, you should set the objective to zero and reoptimize.",
    ),
    (
        MOI.OBJECTIVE_LIMIT,
        "Optimal objective for model was proven to be worse than the value specified in the Cutoff parameter. No solution information is available.",
    ),
    (
        MOI.ITERATION_LIMIT,
        "Optimization terminated because the total number of simplex iterations performed exceeded the value specified in the IterationLimit parameter, or because the total number of barrier iterations exceeded the value specified in the BarIterLimit parameter.",
    ),
    (
        MOI.NODE_LIMIT,
        "Optimization terminated because the total number of branch-and-cut nodes explored exceeded the value specified in the NodeLimit parameter.",
    ),
    (
        MOI.TIME_LIMIT,
        "Optimization terminated because the time expended exceeded the value specified in the TimeLimit parameter.",
    ),
    (
        MOI.SOLUTION_LIMIT,
        "Optimization terminated because the number of solutions found reached the value specified in the SolutionLimit parameter.",
    ),
    (MOI.INTERRUPTED, "Optimization was terminated by the user."),
    (
        MOI.NUMERICAL_ERROR,
        "Optimization was terminated due to unrecoverable numerical difficulties.",
    ),
    (
        MOI.LOCALLY_SOLVED,
        "Unable to satisfy optimality tolerances; a sub-optimal solution is available.",
    ),
    (
        MOI.OTHER_ERROR,
        "An asynchronous optimization call was made, but the associated optimization run is not yet complete.",
    ),
    (
        MOI.OBJECTIVE_LIMIT,
        "User specified an objective limit (a bound on either the best objective or the best bound), and that limit has been reached.",
    ),
]

function _raw_status(model::Optimizer)
    if haskey(_ERROR_TO_STATUS, model.ret_GRBoptimize)
        return _ERROR_TO_STATUS[model.ret_GRBoptimize]
    end
    valueP = Ref{Cint}()
    ret = GRBgetintattr(model, "Status", valueP)
    _check_ret(model, ret)
    @assert 1 <= valueP[] <= 15
    return _RAW_STATUS_STRINGS[valueP[]]
end

function MOI.get(model::Optimizer, attr::MOI.RawStatusString)
    _throw_if_optimize_in_progress(model, attr)
    return _raw_status(model)[2]
end

function MOI.get(model::Optimizer, attr::MOI.TerminationStatus)
    _throw_if_optimize_in_progress(model, attr)
    return _raw_status(model)[1]
end

function _get_intattr(model, key)
    p = Ref{Cint}()
    ret = GRBgetintattr(model, key, p)
    _check_ret(model, ret)
    return p[]
end

function MOI.get(model::Optimizer, attr::MOI.PrimalStatus)
    _throw_if_optimize_in_progress(model, attr)
    term = MOI.get(model, MOI.TerminationStatus())
    if term == MOI.DUAL_INFEASIBLE || term == MOI.INFEASIBLE_OR_UNBOUNDED
        if attr.result_index != 1
            return MOI.NO_SOLUTION
        elseif _has_primal_ray(model)
            return MOI.INFEASIBILITY_CERTIFICATE
        else
            return MOI.NO_SOLUTION
        end
    end
    # Check SolCount explicitly instead of ResultCount to avoid returning 1 when
    # there is a certiticate.
    valueP = Ref{Cint}()
    ret = GRBgetintattr(model, "SolCount", valueP)
    _check_ret(model, ret)
    return 1 <= attr.result_index <= valueP[] ? MOI.FEASIBLE_POINT :
           MOI.NO_SOLUTION
end

function _has_dual_ray(model::Optimizer)
    # Note: for performance reasons, we try to get 1 element because for
    # some versions of Gurobi, we cannot query 0 elements without error.
    valueP = Ref{Cdouble}()
    ret = GRBgetdblattrarray(model, "FarkasDual", 0, 1, valueP)
    return ret == 0
end

_is_mip(model::Optimizer) = _get_intattr(model, "IsMIP") != 0
_is_qcp(model::Optimizer) = _get_intattr(model, "IsQCP") != 0

function MOI.get(model::Optimizer, attr::MOI.DualStatus)
    _throw_if_optimize_in_progress(model, attr)
    if attr.result_index != 1
        return MOI.NO_SOLUTION
    elseif _is_mip(model)
        return MOI.NO_SOLUTION
    elseif _is_qcp(model) &&
           MOI.get(model, MOI.RawOptimizerAttribute("QCPDual")) != 1
        return MOI.NO_SOLUTION
    end
    term = MOI.get(model, MOI.TerminationStatus())
    if term == MOI.INFEASIBLE || term == MOI.INFEASIBLE_OR_UNBOUNDED
        if _has_dual_ray(model)
            return MOI.INFEASIBILITY_CERTIFICATE
        else
            return MOI.NO_SOLUTION
        end
    end
    # Check SolCount explicitly instead of ResultCount to avoid returning 1 when
    # there is a certiticate.
    valueP = Ref{Cint}()
    ret = GRBgetintattr(model, "SolCount", valueP)
    _check_ret(model, ret)
    return valueP[] > 0 ? MOI.FEASIBLE_POINT : MOI.NO_SOLUTION
end

function _has_primal_ray(model::Optimizer)
    # Note: for performance reasons, we try to get 1 element because for
    # some versions of Gurobi, we cannot query 0 elements without error.
    valueP = Ref{Cdouble}()
    ret = GRBgetdblattrarray(model, "UnbdRay", 0, 1, valueP)
    return ret == 0
end

function MOI.get(
    model::Optimizer,
    attr::MOI.VariablePrimal,
    x::MOI.VariableIndex,
)
    _throw_if_optimize_in_progress(model, attr)
    MOI.check_result_index_bounds(model, attr)
    key = model.has_unbounded_ray ? "UnbdRay" : "X"
    if attr.result_index > 1
        MOI.set(
            model,
            MOI.RawOptimizerAttribute("SolutionNumber"),
            attr.result_index - 1,
        )
        key = "Xn"
    end
    col = Cint(column(model, x) - 1)
    valueP = Ref{Cdouble}()
    ret = GRBgetdblattrelement(model, key, col, valueP)
    _check_ret(model, ret)
    return valueP[]
end

function MOI.get(
    model::Optimizer,
    attr::MOI.ConstraintPrimal,
    c::MOI.ConstraintIndex{MOI.VariableIndex,<:Any},
)
    _throw_if_optimize_in_progress(model, attr)
    MOI.check_result_index_bounds(model, attr)
    return MOI.get(model, MOI.VariablePrimal(), MOI.VariableIndex(c.value))
end

function MOI.get(
    model::Optimizer,
    attr::MOI.ConstraintPrimal,
    c::MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64},<:Any},
)
    _throw_if_optimize_in_progress(model, attr)
    MOI.check_result_index_bounds(model, attr)
    row = _info(model, c).row
    _update_if_necessary(model)
    rhs = Ref{Cdouble}()
    ret = GRBgetdblattrelement(model, "RHS", Cint(row - 1), rhs)
    _check_ret(model, ret)
    slack = Ref{Cdouble}()
    ret = GRBgetdblattrelement(model, "Slack", Cint(row - 1), slack)
    _check_ret(model, ret)
    return rhs[] - slack[]
end

function MOI.get(
    model::Optimizer,
    attr::MOI.ConstraintPrimal,
    c::MOI.ConstraintIndex{MOI.ScalarQuadraticFunction{Float64},<:Any},
)
    _throw_if_optimize_in_progress(model, attr)
    MOI.check_result_index_bounds(model, attr)
    row = Cint(_info(model, c).row - 1)
    _update_if_necessary(model)
    rhs = Ref{Cdouble}()
    ret = GRBgetdblattrelement(model, "QCRHS", row, rhs)
    _check_ret(model, ret)
    slack = Ref{Cdouble}()
    ret = GRBgetdblattrelement(model, "QCSlack", row, slack)
    _check_ret(model, ret)
    return rhs[] - slack[]
end

function _dual_multiplier(model::Optimizer)
    return MOI.get(model, MOI.ObjectiveSense()) == MOI.MIN_SENSE ? 1.0 : -1.0
end

"""
    _farkas_variable_dual(model::Optimizer, col::Cint)

Return a Farkas dual associated with the variable bounds of `col`.

Gurobi computes the Farkas dual as:

    ā * x = λ' * A * x <= λ' * b = -β + sum(āᵢ * Uᵢ | āᵢ < 0) + sum(āᵢ * Lᵢ | āᵢ > 0)

The Farkas dual of the variable is ā, and it applies to the upper bound if ā < 0,
and it applies to the lower bound if ā > 0.
"""
function _farkas_variable_dual(model::Optimizer, col::Cint)
    numnzP = Ref{Cint}()
    ret = GRBgetvars(model, numnzP, C_NULL, C_NULL, C_NULL, col, 1)
    _check_ret(model, ret)
    vbeg = Vector{Cint}(undef, 2)
    vind = Vector{Cint}(undef, numnzP[])
    vval = Vector{Cdouble}(undef, numnzP[])
    ret = GRBgetvars(model, numnzP, vbeg, vind, vval, col, 1)
    _check_ret(model, ret)
    λ = Vector{Cdouble}(undef, numnzP[])
    ret = GRBgetdblattrlist(model, "FarkasDual", length(vind), vind, λ)
    _check_ret(model, ret)
    return λ' * vval
end

function MOI.get(
    model::Optimizer,
    attr::MOI.ConstraintDual,
    c::MOI.ConstraintIndex{MOI.VariableIndex,MOI.LessThan{Float64}},
)
    _throw_if_optimize_in_progress(model, attr)
    MOI.check_result_index_bounds(model, attr)
    col = Cint(column(model, c) - 1)
    if model.has_infeasibility_cert
        dual = _farkas_variable_dual(model, col)
        return min(dual, 0.0)
    end
    reduced_cost = Ref{Cdouble}()
    ret = GRBgetdblattrelement(model, "RC", col, reduced_cost)
    _check_ret(model, ret)
    sense = MOI.get(model, MOI.ObjectiveSense())
    # The following is a heuristic for determining whether the reduced cost
    # applies to the lower or upper bound. It can be wrong by at most
    # `FeasibilityTol`.
    if sense == MOI.MIN_SENSE && reduced_cost[] < 0
        # If minimizing, the reduced cost must be negative (ignoring
        # tolerances).
        return reduced_cost[]
    elseif sense == MOI.MAX_SENSE && reduced_cost[] > 0
        # If maximizing, the reduced cost must be positive (ignoring
        # tolerances). However, because of the MOI dual convention,
        # we return a negative value.
        return -reduced_cost[]
    else
        # The reduced cost, if non-zero, must related to the lower bound.
        return 0.0
    end
end

function MOI.get(
    model::Optimizer,
    attr::MOI.ConstraintDual,
    c::MOI.ConstraintIndex{MOI.VariableIndex,MOI.GreaterThan{Float64}},
)
    _throw_if_optimize_in_progress(model, attr)
    MOI.check_result_index_bounds(model, attr)
    col = Cint(column(model, c) - 1)
    if model.has_infeasibility_cert
        dual = _farkas_variable_dual(model, col)
        return max(dual, 0.0)
    end
    reduced_cost = Ref{Cdouble}()
    ret = GRBgetdblattrelement(model, "RC", col, reduced_cost)
    _check_ret(model, ret)
    sense = MOI.get(model, MOI.ObjectiveSense())
    # The following is a heuristic for determining whether the reduced cost
    # applies to the lower or upper bound. It can be wrong by at most
    # `FeasibilityTol`.
    if sense == MOI.MIN_SENSE && reduced_cost[] > 0
        # If minimizing, the reduced cost must be positive (ignoring
        # tolerances).
        return reduced_cost[]
    elseif sense == MOI.MAX_SENSE && reduced_cost[] < 0
        # If maximizing, the reduced cost must be negative (ignoring
        # tolerances). However, because of the MOI dual convention,
        # we return a positive value.
        return -reduced_cost[]
    else
        # The reduced cost, if non-zero, must related to the lower bound.
        return 0.0
    end
end

function MOI.get(
    model::Optimizer,
    attr::MOI.ConstraintDual,
    c::MOI.ConstraintIndex{MOI.VariableIndex,MOI.EqualTo{Float64}},
)
    _throw_if_optimize_in_progress(model, attr)
    MOI.check_result_index_bounds(model, attr)
    col = Cint(column(model, c) - 1)
    if model.has_infeasibility_cert
        return _farkas_variable_dual(model, col)
    end
    reduced_cost = Ref{Cdouble}()
    ret = GRBgetdblattrelement(model, "RC", col, reduced_cost)
    _check_ret(model, ret)
    return _dual_multiplier(model) * reduced_cost[]
end

function MOI.get(
    model::Optimizer,
    attr::MOI.ConstraintDual,
    c::MOI.ConstraintIndex{MOI.VariableIndex,MOI.Interval{Float64}},
)
    _throw_if_optimize_in_progress(model, attr)
    MOI.check_result_index_bounds(model, attr)
    col = Cint(column(model, c) - 1)
    if model.has_infeasibility_cert
        return _farkas_variable_dual(model, col)
    end
    reduced_cost = Ref{Cdouble}()
    ret = GRBgetdblattrelement(model, "RC", col, reduced_cost)
    _check_ret(model, ret)
    return _dual_multiplier(model) * reduced_cost[]
end

function MOI.get(
    model::Optimizer,
    attr::MOI.ConstraintDual,
    c::MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64},<:Any},
)
    _throw_if_optimize_in_progress(model, attr)
    MOI.check_result_index_bounds(model, attr)
    valueP = Ref{Cdouble}()
    row = Cint(_info(model, c).row - 1)
    if model.has_infeasibility_cert
        ret = GRBgetdblattrelement(model, "FarkasDual", row, valueP)
        _check_ret(model, ret)
        return -valueP[]
    end
    ret = GRBgetdblattrelement(model, "Pi", row, valueP)
    _check_ret(model, ret)
    return _dual_multiplier(model) * valueP[]
end

function MOI.get(
    model::Optimizer,
    attr::MOI.ConstraintDual,
    c::MOI.ConstraintIndex{MOI.ScalarQuadraticFunction{Float64},<:Any},
)
    _throw_if_optimize_in_progress(model, attr)
    MOI.check_result_index_bounds(model, attr)
    if model.has_infeasibility_cert
        error("Infeasibility certificate not available for $(c).")
    end
    valueP = Ref{Cdouble}()
    ret = GRBgetdblattrelement(
        model,
        "QCPi",
        Cint(_info(model, c).row - 1),
        valueP,
    )
    _check_ret(model, ret)
    return _dual_multiplier(model) * valueP[]
end

function MOI.get(model::Optimizer, attr::MOI.ObjectiveValue)
    _throw_if_optimize_in_progress(model, attr)
    MOI.check_result_index_bounds(model, attr)
    if attr.result_index > 1
        MOI.set(
            model,
            MOI.RawOptimizerAttribute("SolutionNumber"),
            attr.result_index - 1,
        )
    end
    valueP = Ref{Cdouble}()
    key = attr.result_index == 1 ? "ObjVal" : "PoolObjVal"
    ret = GRBgetdblattr(model, key, valueP)
    _check_ret(model, ret)
    return valueP[]
end

function MOI.get(model::Optimizer, attr::MOI.ObjectiveBound)
    _throw_if_optimize_in_progress(model, attr)
    valueP = Ref{Cdouble}()
    ret = GRBgetdblattr(model, "ObjBound", valueP)
    _check_ret(model, ret)
    return valueP[]
end

function MOI.get(model::Optimizer, attr::MOI.SolveTimeSec)
    _throw_if_optimize_in_progress(model, attr)
    valueP = Ref{Cdouble}()
    ret = GRBgetdblattr(model, "RunTime", valueP)
    _check_ret(model, ret)
    return valueP[]
end

function MOI.get(model::Optimizer, attr::MOI.SimplexIterations)
    _throw_if_optimize_in_progress(model, attr)
    valueP = Ref{Cdouble}()
    ret = GRBgetdblattr(model, "IterCount", valueP)
    _check_ret(model, ret)
    return valueP[]
end

function MOI.get(model::Optimizer, attr::MOI.BarrierIterations)
    _throw_if_optimize_in_progress(model, attr)
    valueP = Ref{Cint}()
    ret = GRBgetintattr(model, "BarIterCount", valueP)
    _check_ret(model, ret)
    return Int(valueP[])
end

function MOI.get(model::Optimizer, attr::MOI.NodeCount)
    _throw_if_optimize_in_progress(model, attr)
    valueP = Ref{Cdouble}()
    ret = GRBgetdblattr(model, "NodeCount", valueP)
    _check_ret(model, ret)
    return valueP[]
end

function MOI.get(model::Optimizer, attr::MOI.RelativeGap)
    _throw_if_optimize_in_progress(model, attr)
    valueP = Ref{Cdouble}()
    ret = GRBgetdblattr(model, "MIPGap", valueP)
    _check_ret(model, ret)
    return valueP[]
end

function MOI.get(model::Optimizer, attr::MOI.DualObjectiveValue)
    _throw_if_optimize_in_progress(model, attr)
    MOI.check_result_index_bounds(model, attr)
    valueP = Ref{Cdouble}()
    ret = GRBgetdblattr(model, "ObjBound", valueP)
    _check_ret(model, ret)
    return valueP[]
end

function MOI.get(model::Optimizer, attr::MOI.ResultCount)
    _throw_if_optimize_in_progress(model, attr)
    if model.has_infeasibility_cert || model.has_unbounded_ray
        return 1
    end
    valueP = Ref{Cint}()
    ret = GRBgetintattr(model, "SolCount", valueP)
    _check_ret(model, ret)
    return Int(valueP[])
end

function MOI.get(model::Optimizer, ::MOI.Silent)
    return model.silent
end

function MOI.set(model::Optimizer, ::MOI.Silent, flag::Bool)
    model.silent = flag
    output_flag = flag ? 0 : 1
    MOI.set(model, MOI.RawOptimizerAttribute("OutputFlag"), output_flag)
    return
end

function MOI.get(model::Optimizer, ::MOI.NumberOfThreads)
    x = MOI.get(model, MOI.RawOptimizerAttribute("Threads"))
    # Instead of default `0`, return `nothing`
    return x == 0 ? nothing : x
end

function MOI.set(model::Optimizer, ::MOI.NumberOfThreads, x::Int)
    MOI.set(model, MOI.RawOptimizerAttribute("Threads"), x)
    return
end

function MOI.set(model::Optimizer, ::MOI.NumberOfThreads, ::Nothing)
    MOI.set(model, MOI.RawOptimizerAttribute("Threads"), 0)
    return
end

function MOI.get(model::Optimizer, ::MOI.Name)
    _update_if_necessary(model)
    valueP = Ref{Ptr{Cchar}}()
    ret = GRBgetstrattr(model, "ModelName", valueP)
    _check_ret(model, ret)
    return unsafe_string(valueP[])
end

function MOI.set(model::Optimizer, ::MOI.Name, name::String)
    ret = GRBsetstrattr(model, "ModelName", name)
    _check_ret(model, ret)
    _require_update(model)
    return
end

MOI.get(model::Optimizer, ::MOI.NumberOfVariables) = length(model.variable_info)
function MOI.get(model::Optimizer, ::MOI.ListOfVariableIndices)
    return sort!(collect(keys(model.variable_info)), by = x -> x.value)
end

MOI.get(model::Optimizer, ::MOI.RawSolver) = model

function MOI.set(
    model::Optimizer,
    ::MOI.VariablePrimalStart,
    x::MOI.VariableIndex,
    value::Union{Nothing,Float64},
)
    info = _info(model, x)
    info.start = value
    grb_value = value !== nothing ? value : GRB_UNDEFINED
    ret = GRBsetdblattrelement(model, "Start", Cint(info.column - 1), grb_value)
    _check_ret(model, ret)
    _require_update(model)
    return
end

function MOI.get(
    model::Optimizer,
    ::MOI.VariablePrimalStart,
    x::MOI.VariableIndex,
)
    return _info(model, x).start
end

function MOI.supports(
    ::Gurobi.Optimizer,
    ::MOI.VariablePrimalStart,
    ::Type{MOI.VariableIndex},
)
    return true
end

function MOI.get(model::Optimizer, ::MOI.NumberOfConstraints{F,S}) where {F,S}
    # TODO: this could be more efficient.
    return length(MOI.get(model, MOI.ListOfConstraintIndices{F,S}()))
end

_bound_enums(::Type{<:MOI.LessThan}) = (_LESS_THAN, _LESS_AND_GREATER_THAN)
function _bound_enums(::Type{<:MOI.GreaterThan})
    return (_GREATER_THAN, _LESS_AND_GREATER_THAN)
end
_bound_enums(::Type{<:MOI.Interval}) = (_INTERVAL,)
_bound_enums(::Type{<:MOI.EqualTo}) = (_EQUAL_TO,)
_bound_enums(::Any) = (nothing,)

_type_enums(::Type{MOI.ZeroOne}) = (GRB_BINARY,)
_type_enums(::Type{MOI.Integer}) = (GRB_INTEGER,)
_type_enums(::Type{<:MOI.Semicontinuous}) = (GRB_SEMICONT,)
_type_enums(::Type{<:MOI.Semiinteger}) = (GRB_SEMIINT,)
_type_enums(::Any) = (nothing,)

function MOI.get(
    model::Optimizer,
    ::MOI.ListOfConstraintIndices{MOI.VariableIndex,S},
) where {S}
    indices = MOI.ConstraintIndex{MOI.VariableIndex,S}[]
    for (key, info) in model.variable_info
        if info.bound in _bound_enums(S) || info.type in _type_enums(S)
            push!(indices, MOI.ConstraintIndex{MOI.VariableIndex,S}(key.value))
        end
    end
    return sort!(indices, by = x -> x.value)
end

function MOI.get(
    model::Optimizer,
    ::MOI.ListOfConstraintIndices{MOI.ScalarAffineFunction{Float64},S},
) where {S}
    indices = MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64},S}[]
    for (key, info) in model.affine_constraint_info
        if typeof(info.set) == S
            push!(
                indices,
                MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64},S}(key),
            )
        end
    end
    return sort!(indices, by = x -> x.value)
end

function MOI.get(
    model::Optimizer,
    ::MOI.ListOfConstraintIndices{MOI.ScalarQuadraticFunction{Float64},S},
) where {S}
    indices = MOI.ConstraintIndex{MOI.ScalarQuadraticFunction{Float64},S}[]
    for (key, info) in model.quadratic_constraint_info
        if typeof(info.set) == S
            push!(
                indices,
                MOI.ConstraintIndex{MOI.ScalarQuadraticFunction{Float64},S}(
                    key,
                ),
            )
        end
    end
    return sort!(indices, by = x -> x.value)
end

function MOI.get(
    model::Optimizer,
    ::MOI.ListOfConstraintIndices{MOI.VectorOfVariables,S},
) where {S<:Union{<:MOI.SOS1,<:MOI.SOS2}}
    indices = MOI.ConstraintIndex{MOI.VectorOfVariables,S}[]
    for (key, info) in model.sos_constraint_info
        if typeof(info.set) == S
            push!(indices, MOI.ConstraintIndex{MOI.VectorOfVariables,S}(key))
        end
    end
    return sort!(indices, by = x -> x.value)
end

function MOI.get(
    model::Optimizer,
    ::MOI.ListOfConstraintIndices{MOI.VectorOfVariables,MOI.SecondOrderCone},
)
    indices = MOI.ConstraintIndex{MOI.VectorOfVariables,MOI.SecondOrderCone}[
        MOI.ConstraintIndex{MOI.VectorOfVariables,MOI.SecondOrderCone}(key)
        for (key, info) in model.quadratic_constraint_info if
        typeof(info.set) == MOI.SecondOrderCone
    ]
    return sort!(indices, by = x -> x.value)
end

function MOI.get(model::Optimizer, ::MOI.ListOfConstraintTypesPresent)
    constraints = Set{Tuple{Type,Type}}()
    for info in values(model.variable_info)
        if info.bound == _NONE
        elseif info.bound == _LESS_THAN
            push!(constraints, (MOI.VariableIndex, MOI.LessThan{Float64}))
        elseif info.bound == _GREATER_THAN
            push!(constraints, (MOI.VariableIndex, MOI.GreaterThan{Float64}))
        elseif info.bound == _LESS_AND_GREATER_THAN
            push!(constraints, (MOI.VariableIndex, MOI.LessThan{Float64}))
            push!(constraints, (MOI.VariableIndex, MOI.GreaterThan{Float64}))
        elseif info.bound == _EQUAL_TO
            push!(constraints, (MOI.VariableIndex, MOI.EqualTo{Float64}))
        elseif info.bound == _INTERVAL
            push!(constraints, (MOI.VariableIndex, MOI.Interval{Float64}))
        end
        if info.type == GRB_CONTINUOUS
        elseif info.type == GRB_BINARY
            push!(constraints, (MOI.VariableIndex, MOI.ZeroOne))
        elseif info.type == GRB_INTEGER
            push!(constraints, (MOI.VariableIndex, MOI.Integer))
        elseif info.type == GRB_SEMICONT
            push!(constraints, (MOI.VariableIndex, MOI.Semicontinuous{Float64}))
        elseif info.type == GRB_SEMIINT
            push!(constraints, (MOI.VariableIndex, MOI.Semiinteger{Float64}))
        end
    end
    for info in values(model.affine_constraint_info)
        push!(
            constraints,
            (MOI.ScalarAffineFunction{Float64}, typeof(info.set)),
        )
    end
    for info in values(model.quadratic_constraint_info)
        if typeof(info.set) == MOI.SecondOrderCone
            push!(constraints, (MOI.VectorOfVariables, MOI.SecondOrderCone))
        else
            push!(
                constraints,
                (MOI.ScalarQuadraticFunction{Float64}, typeof(info.set)),
            )
        end
    end
    for info in values(model.sos_constraint_info)
        push!(constraints, (MOI.VectorOfVariables, typeof(info.set)))
    end
    return collect(constraints)
end

function MOI.get(model::Optimizer, ::MOI.ObjectiveFunctionType)
    if model.objective_type == _SINGLE_VARIABLE
        return MOI.VariableIndex
    elseif model.objective_type == _SCALAR_AFFINE
        return MOI.ScalarAffineFunction{Float64}
    else
        @assert model.objective_type == _SCALAR_QUADRATIC
        return MOI.ScalarQuadraticFunction{Float64}
    end
end

function MOI.modify(
    model::Optimizer,
    c::MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64},<:Any},
    chg::MOI.ScalarCoefficientChange{Float64},
)
    ret = GRBchgcoeffs(
        model,
        1,
        Ref{Cint}(_info(model, c).row - 1),
        Ref{Cint}(column(model, chg.variable) - 1),
        Ref{Cdouble}(chg.new_coefficient),
    )
    _check_ret(model, ret)
    _require_update(model)
    return
end

function MOI.modify(
    model::Optimizer,
    ::MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}},
    chg::MOI.ScalarCoefficientChange{Float64},
)
    ret = GRBsetdblattrelement(
        model,
        "Obj",
        Cint(column(model, chg.variable) - 1),
        chg.new_coefficient,
    )
    _check_ret(model, ret)
    model.is_objective_set = true
    _require_update(model)
    return
end

"""
    _replace_with_matching_sparsity!(
        model::Optimizer,
        previous::MOI.ScalarAffineFunction,
        replacement::MOI.ScalarAffineFunction, row::Int
    )

Internal function, not intended for external use.

Change the linear constraint function at index `row` in `model` from
`previous` to `replacement`. This function assumes that `previous` and
`replacement` have exactly the same sparsity pattern w.r.t. which variables
they include and that both constraint functions are in canonical form (as
returned by `MOIU.canonical()`. Neither assumption is checked within the body
of this function.
"""
function _replace_with_matching_sparsity!(
    model::Optimizer,
    previous::MOI.ScalarAffineFunction,
    replacement::MOI.ScalarAffineFunction,
    row::Int,
)
    rows = fill(Cint(row - 1), length(replacement.terms))
    cols = Cint[column(model, t.variable) - 1 for t in replacement.terms]
    coefs = MOI.coefficient.(replacement.terms)
    ret = GRBchgcoeffs(model, length(cols), rows, cols, coefs)
    _check_ret(model, ret)
    return
end

"""
    _replace_with_different_sparsity!(
        model::Optimizer,
        previous::MOI.ScalarAffineFunction,
        replacement::MOI.ScalarAffineFunction, row::Int
    )

Internal function, not intended for external use.

    Change the linear constraint function at index `row` in `model` from
`previous` to `replacement`. This function assumes that `previous` and
`replacement` may have different sparsity patterns.

This function (and `_replace_with_matching_sparsity!` above) are necessary
because in order to fully replace a linear constraint, we have to zero out the
current matrix coefficients and then set the new matrix coefficients. When the
sparsity patterns match, the zeroing-out step can be skipped.
"""
function _replace_with_different_sparsity!(
    model::Optimizer,
    previous::MOI.ScalarAffineFunction,
    replacement::MOI.ScalarAffineFunction,
    row::Int,
)
    # First, zero out the old constraint function terms.
    rows = fill(Cint(row - 1), length(previous.terms))
    cols = Cint[column(model, t.variable) - 1 for t in previous.terms]
    coefs = fill(0.0, length(previous.terms))
    ret = GRBchgcoeffs(model, length(cols), rows, cols, coefs)
    _check_ret(model, ret)
    # Next, set the new constraint function terms.
    rows = fill(Cint(row - 1), length(replacement.terms))
    cols = Cint[column(model, t.variable) - 1 for t in replacement.terms]
    coefs = MOI.coefficient.(replacement.terms)
    ret = GRBchgcoeffs(model, length(cols), rows, cols, coefs)
    _check_ret(model, ret)
    return
end

"""
    _matching_sparsity_pattern(
        f1::MOI.ScalarAffineFunction{Float64},
        f2::MOI.ScalarAffineFunction{Float64}
    )

Internal function, not intended for external use.

Determines whether functions `f1` and `f2` have the same sparsity pattern
w.r.t. their constraint columns. Assumes both functions are already in
canonical form.
"""
function _matching_sparsity_pattern(
    f1::MOI.ScalarAffineFunction{Float64},
    f2::MOI.ScalarAffineFunction{Float64},
)
    if axes(f1.terms) != axes(f2.terms)
        return false
    end
    for (f1_term, f2_term) in zip(f1.terms, f2.terms)
        if MOI.term_indices(f1_term) != MOI.term_indices(f2_term)
            return false
        end
    end
    return true
end

function MOI.set(
    model::Optimizer,
    ::MOI.ConstraintFunction,
    c::MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64},<:_SCALAR_SETS},
    f::MOI.ScalarAffineFunction{Float64},
)
    previous = MOI.get(model, MOI.ConstraintFunction(), c)
    MOI.Utilities.canonicalize!(previous)
    replacement = MOI.Utilities.canonical(f)
    _update_if_necessary(model)
    # If the previous and replacement constraint functions have exactly
    # the same sparsity pattern, then we can take a faster path by just
    # passing the replacement terms to the model. But if their sparsity
    # patterns differ, then we need to first zero out the previous terms
    # and then set the replacement terms.
    row = _info(model, c).row
    if _matching_sparsity_pattern(previous, replacement)
        _replace_with_matching_sparsity!(model, previous, replacement, row)
    else
        _replace_with_different_sparsity!(model, previous, replacement, row)
    end
    current_rhs = Ref{Cdouble}()
    ret = GRBgetdblattrelement(model, "RHS", Cint(row - 1), current_rhs)
    _check_ret(model, ret)
    new_rhs = current_rhs[] - (replacement.constant - previous.constant)
    ret = GRBsetdblattrelement(model, "RHS", Cint(row - 1), new_rhs)
    _check_ret(model, ret)
    _require_update(model)
    return
end

function MOI.get(
    model::Optimizer,
    ::MOI.ConstraintBasisStatus,
    c::MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64},S},
) where {S<:_SCALAR_SETS}
    row = _info(model, c).row
    _update_if_necessary(model)
    valueP = Ref{Cint}()
    ret = GRBgetintattrelement(model, "CBasis", Cint(row - 1), valueP)
    _check_ret(model, ret)
    if valueP[] == 0
        return MOI.BASIC
    elseif valueP[] == -1
        return MOI.NONBASIC
    else
        error("CBasis value of $(valueP[]) isn't defined.")
    end
end

function MOI.get(
    model::Optimizer,
    ::MOI.VariableBasisStatus,
    x::MOI.VariableIndex,
)
    _update_if_necessary(model)
    valueP = Ref{Cint}()
    col = Cint(column(model, x) - 1)
    ret = GRBgetintattrelement(model, "VBasis", col, valueP)
    _check_ret(model, ret)
    if valueP[] == 0
        return MOI.BASIC
    elseif valueP[] == -1
        return MOI.NONBASIC_AT_LOWER
    elseif valueP[] == -2
        return MOI.NONBASIC_AT_UPPER
    else
        @assert valueP[] == -3
        return MOI.SUPER_BASIC
    end
end

function MOI.set(
    model::Optimizer,
    ::MOI.ConstraintBasisStatus,
    c::MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64},S},
    bs::MOI.BasisStatusCode,
) where {S<:_SCALAR_SETS}
    row = _info(model, c).row
    valueP::Cint = (
        if bs == MOI.BASIC
            Cint(-0)
        elseif bs == MOI.NONBASIC
            Cint(-1)
        else
            error("Unsupported constraint basis value $bs.")
        end
    )
    ret = GRBsetintattrelement(model, "CBasis", Cint(row - 1), valueP)
    _check_ret(model, ret)
    _require_update(model)
    return
end

function MOI.set(
    model::Optimizer,
    ::MOI.VariableBasisStatus,
    x::MOI.VariableIndex,
    bs::MOI.BasisStatusCode,
)
    valueP = if bs == MOI.BASIC
        Cint(0)
    elseif bs == MOI.NONBASIC_AT_LOWER
        Cint(-1)
    elseif bs == MOI.NONBASIC_AT_UPPER
        Cint(-2)
    else
        @assert bs == MOI.SUPER_BASIC
        Cint(-3)
    end
    col = Cint(column(model, x) - 1)
    ret = GRBsetintattrelement(model, "VBasis", col, valueP)
    _check_ret(model, ret)
    _require_update(model)
    return
end

function MOI.compute_conflict!(model::Optimizer)
    ret = GRBcomputeIIS(model)
    _check_ret(model, ret)
    model.conflict = ret
    return
end

function _ensure_conflict_computed(model::Optimizer)
    if model.conflict == -1
        error(
            "Cannot access conflict status. Call " *
            "`Gurobi.compute_conflict(model)` first. In case the model " *
            "is modified, the computed conflict will not be purged.",
        )
    end
end

function _is_feasible(model::Optimizer)
    return model.conflict == Gurobi.GRB_INFEASIBLE
end

"""
    ConflictStatus()

Return the raw status from Gurobi indicating the status of the last
computed conflict. It returns an integer:

* `-1` if `compute_conflict!` has not yet been called
* `0` if it found a conflict
* other values are defined in [Gurobi's documentation](https://www.gurobi.com/documentation/current/refman/error_codes.html)
"""
struct ConflictStatus <: MOI.AbstractModelAttribute end

function MOI.get(model::Optimizer, ::ConflictStatus)
    return model.conflict
end

function MOI.get(model::Optimizer, ::MOI.ConflictStatus)
    status = MOI.get(model, ConflictStatus())
    if status == -1
        return MOI.COMPUTE_CONFLICT_NOT_CALLED
    elseif status == 0
        return MOI.CONFLICT_FOUND
    elseif status == Gurobi.IIS_NOT_INFEASIBLE
        return MOI.NO_CONFLICT_EXISTS
    else
        return MOI.NO_CONFLICT_FOUND
    end
end

function MOI.get(
    model::Optimizer,
    ::MOI.ConstraintConflictStatus,
    index::MOI.ConstraintIndex{MOI.VariableIndex,<:MOI.LessThan},
)
    _ensure_conflict_computed(model)
    if _is_feasible(model)
        return MOI.NOT_IN_CONFLICT
    end
    p = Ref{Cint}()
    ret =
        GRBgetintattrelement(model, "IISUB", Cint(column(model, index) - 1), p)
    _check_ret(model, ret)
    return p[] > 0 ? MOI.IN_CONFLICT : MOI.NOT_IN_CONFLICT
end

function MOI.get(
    model::Optimizer,
    ::MOI.ConstraintConflictStatus,
    index::MOI.ConstraintIndex{MOI.VariableIndex,<:MOI.GreaterThan},
)
    _ensure_conflict_computed(model)
    if _is_feasible(model)
        return MOI.NOT_IN_CONFLICT
    end
    p = Ref{Cint}()
    ret =
        GRBgetintattrelement(model, "IISLB", Cint(column(model, index) - 1), p)
    _check_ret(model, ret)
    return p[] > 0 ? MOI.IN_CONFLICT : MOI.NOT_IN_CONFLICT
end

function MOI.get(
    model::Optimizer,
    ::MOI.ConstraintConflictStatus,
    index::MOI.ConstraintIndex{
        MOI.VariableIndex,
        <:Union{MOI.EqualTo,MOI.Interval},
    },
)
    _ensure_conflict_computed(model)
    if _is_feasible(model)
        return MOI.NOT_IN_CONFLICT
    end
    p = Ref{Cint}()
    ret =
        GRBgetintattrelement(model, "IISLB", Cint(column(model, index) - 1), p)
    _check_ret(model, ret)
    if p[] > 0
        return MOI.IN_CONFLICT
    end
    ret =
        GRBgetintattrelement(model, "IISUB", Cint(column(model, index) - 1), p)
    _check_ret(model, ret)
    if p[] > 0
        return MOI.IN_CONFLICT
    end
    return MOI.NOT_IN_CONFLICT
end

function MOI.get(
    model::Optimizer,
    ::MOI.ConstraintConflictStatus,
    index::MOI.ConstraintIndex{
        MOI.ScalarAffineFunction{Float64},
        <:_SUPPORTED_SCALAR_SETS,
    },
)
    _ensure_conflict_computed(model)
    if _is_feasible(model)
        return MOI.NOT_IN_CONFLICT
    end
    p = Ref{Cint}()
    ret = GRBgetintattrelement(
        model,
        "IISConstr",
        Cint(_info(model, index).row - 1),
        p,
    )
    _check_ret(model, ret)
    return p[] > 0 ? MOI.IN_CONFLICT : MOI.NOT_IN_CONFLICT
end

function MOI.get(
    model::Optimizer,
    ::MOI.ConstraintConflictStatus,
    index::MOI.ConstraintIndex{
        MOI.ScalarQuadraticFunction{Float64},
        <:Union{MOI.LessThan,MOI.GreaterThan},
    },
)
    _ensure_conflict_computed(model)
    if _is_feasible(model)
        return MOI.NOT_IN_CONFLICT
    end
    p = Ref{Cint}()
    ret = GRBgetintattrelement(
        model,
        "IISQConstr",
        Cint(_info(model, index).row - 1),
        p,
    )
    _check_ret(model, ret)
    return p[] > 0 ? MOI.IN_CONFLICT : MOI.NOT_IN_CONFLICT
end

function MOI.get(
    model::Optimizer,
    ::MOI.ConstraintConflictStatus,
    index::MOI.ConstraintIndex{
        MOI.VariableIndex,
        <:Union{MOI.Integer,MOI.ZeroOne},
    },
)
    _ensure_conflict_computed(model)
    if _is_feasible(model)
        return MOI.NOT_IN_CONFLICT
    end

    # Gurobi doesn't give that information (only linear constraints and bounds,
    # i.e. about the linear relaxation), even though it will report a conflict.
    # Even for binary variables (MOI.ZeroOne), no variable attribute is set by
    # the IIS computation, including the bounds.
    # -> https://www.gurobi.com/documentation/9.1/refman/iislb.html
    # Report that lack of information to the user.
    if MOI.is_valid(model, index)
        return MOI.MAYBE_IN_CONFLICT
    else
        throw(MOI.InvalidIndex(index))
    end
end

###
### Gurobi-specific attributes
###

struct ModelAttribute <: MOI.AbstractModelAttribute
    name::String
end

struct VariableAttribute <: MOI.AbstractVariableAttribute
    name::String
end

struct ConstraintAttribute <: MOI.AbstractConstraintAttribute
    name::String
end

_supported_attrtypes(::ModelAttribute) = (0,)
_supported_attrtypes(::VariableAttribute) = (1,)
_supported_attrtypes(::ConstraintAttribute) = (2, 3, 4, 5)

function _supported(model, attr)
    attrtypeP = Ref{Cint}()
    ret = GRBgetattrinfo(model, attr.name, C_NULL, attrtypeP, C_NULL)
    _check_ret(model, ret)
    return attrtypeP[] in _supported_attrtypes(attr)
end

function MOI.supports(model::Optimizer, attr::ModelAttribute)
    return _supported(model, attr)
end

function MOI.supports(
    model::Optimizer,
    attr::VariableAttribute,
    ::Type{MOI.VariableIndex},
)
    return _supported(model, attr)
end

function MOI.supports(
    model::Optimizer,
    attr::ConstraintAttribute,
    ::Type{<:MOI.ConstraintIndex},
)
    return _supported(model, attr)
end

_check_argument(::String, ::Any, ::Type{S}, ::T) where {S,T<:S} = nothing
function _check_argument(attr::String, T, ::Type{S}, value) where {S}
    return throw(
        ArgumentError(
            "Attribute $(attr) requires $(T) arguments. Provided argument " *
            "was of type $(typeof(value)).",
        ),
    )
end

function _set_attribute(model::Optimizer, attr::ModelAttribute, value)
    name = attr.name
    datatypeP, attrtypeP = Ref{Cint}(), Ref{Cint}()
    ret = GRBgetattrinfo(model, name, datatypeP, attrtypeP, C_NULL)
    _check_ret(model, ret)
    if !(attrtypeP[] in _supported_attrtypes(attr))
        throw(MOI.UnsupportedAttribute(attr))
    end
    ret = if datatypeP[] == 1
        _check_argument(name, Int, Integer, value)
        GRBsetintattr(model, name, value)
    elseif datatypeP[] == 2
        _check_argument(name, Cdouble, Real, value)
        GRBsetdblattr(model, name, value)
    else
        @assert datatypeP[] == 3
        _check_argument(name, String, String, value)
        GRBsetstrattr(model, name, value)
    end
    _check_ret(model, ret)
    return
end

function _get_attribute(model::Optimizer, attr::ModelAttribute)
    name = attr.name
    datatypeP, attrtypeP = Ref{Cint}(), Ref{Cint}()
    ret = GRBgetattrinfo(model, name, datatypeP, attrtypeP, C_NULL)
    _check_ret(model, ret)
    if !(attrtypeP[] in _supported_attrtypes(attr))
        throw(MOI.UnsupportedAttribute(attr))
    end
    if datatypeP[] == 1
        p = Ref{Cint}()
        ret = GRBgetintattr(model, name, p)
        _check_ret(model, ret)
        return Int(p[])
    elseif datatypeP[] == 2
        p = Ref{Cdouble}()
        ret = GRBgetdblattr(model, name, p)
        _check_ret(model, ret)
        return p[]
    else
        @assert datatypeP[] == 3
        valueP = Ref{Ptr{Cchar}}()
        ret = GRBgetstrattr(model, name, valueP)
        _check_ret(model, ret)
        return unsafe_string(valueP[])
    end
    return
end

function _set_attribute(model::Optimizer, attr, element::Cint, value)
    name = attr.name
    datatypeP, attrtypeP = Ref{Cint}(), Ref{Cint}()
    ret = GRBgetattrinfo(model, name, datatypeP, attrtypeP, C_NULL)
    _check_ret(model, ret)
    if !(attrtypeP[] in _supported_attrtypes(attr))
        throw(MOI.UnsupportedAttribute(attr))
    end
    ret = if datatypeP[] == 0
        _check_argument(name, Char, Union{Char,Cchar}, value)
        GRBsetcharattrelement(model, name, element, value)
    elseif datatypeP[] == 1
        _check_argument(name, Int, Integer, value)
        GRBsetintattrelement(model, name, element, value)
    elseif datatypeP[] == 2
        _check_argument(name, Cdouble, Real, value)
        GRBsetdblattrelement(model, name, element, value)
    else
        @assert datatypeP[] == 3
        _check_argument(name, String, String, value)
        GRBsetstrattrelement(model, name, element, value)
    end
    _check_ret(model, ret)
    return
end

function _get_attribute(model::Optimizer, attr, element::Cint)
    name = attr.name
    datatypeP, attrtypeP = Ref{Cint}(), Ref{Cint}()
    ret = GRBgetattrinfo(model, name, datatypeP, attrtypeP, C_NULL)
    _check_ret(model, ret)
    if !(attrtypeP[] in _supported_attrtypes(attr))
        throw(MOI.UnsupportedAttribute(attr))
    end
    if datatypeP[] == 0
        p = Ref{Cchar}()
        ret = GRBgetcharattrelement(model, name, element, p)
        _check_ret(model, ret)
        return Char(p[])
    elseif datatypeP[] == 1
        p = Ref{Cint}()
        ret = GRBgetintattrelement(model, name, element, p)
        _check_ret(model, ret)
        return Int(p[])
    elseif datatypeP[] == 2
        p = Ref{Cdouble}()
        ret = GRBgetdblattrelement(model, name, element, p)
        _check_ret(model, ret)
        return p[]
    else
        @assert datatypeP[] == 3
        valueP = Ref{Ptr{Cchar}}()
        ret = GRBgetstrattrelement(model, name, element, valueP)
        _check_ret(model, ret)
        return unsafe_string(valueP[])
    end
    return
end

"""
    MOI.set(
        model::Optimizer,
        attr::ConstraintAttribute,
        ci::MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64}, <:Any},
        value
    )

Set a constraint attribute.
"""
function MOI.set(
    model::Optimizer,
    attr::ConstraintAttribute,
    ci::MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64}},
    value,
)
    _set_attribute(model, attr, Cint(_info(model, ci).row - 1), value)
    _require_update(model)
    return
end

function MOI.get(
    model::Optimizer,
    attr::ConstraintAttribute,
    ci::MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64}},
)
    _update_if_necessary(model)
    return _get_attribute(model, attr, Cint(_info(model, ci).row - 1))
end

"""
    MOI.set(
        model::Optimizer,
        attr::VariableAttribute,
        vi::MOI.VariableIndex,
        value
    )

Set a variable attribute.
"""
function MOI.set(
    model::Optimizer,
    attr::VariableAttribute,
    vi::MOI.VariableIndex,
    value::T,
) where {T}
    _set_attribute(model, attr, Cint(column(model, vi) - 1), value)
    _require_update(model)
    return
end

function MOI.get(
    model::Optimizer,
    attr::VariableAttribute,
    vi::MOI.VariableIndex,
)
    _update_if_necessary(model)
    return _get_attribute(model, attr, Cint(column(model, vi) - 1))
end

"""
    MOI.set(model::Optimizer, attr::ModelAttribute, value)

Set a model attribute.
"""
function MOI.set(model::Optimizer, attr::ModelAttribute, value)
    _set_attribute(model, attr, value)
    _require_update(model)
    return
end

function MOI.get(model::Optimizer, attr::ModelAttribute)
    _update_if_necessary(model)
    return _get_attribute(model, attr)
end

###
### VectorOfVariables-in-SecondOrderCone
###

function _info(
    model::Optimizer,
    c::MOI.ConstraintIndex{MOI.VectorOfVariables,MOI.SecondOrderCone},
)
    if haskey(model.quadratic_constraint_info, c.value)
        return model.quadratic_constraint_info[c.value]
    end
    return throw(MOI.InvalidIndex(c))
end

function MOI.add_constraint(
    model::Optimizer,
    f::MOI.VectorOfVariables,
    s::MOI.SecondOrderCone,
)
    if length(f.variables) != s.dimension
        error("Dimension of $(s) does not match number of terms in $(f)")
    end

    # SOC is the cone: t ≥ ||x||₂ ≥ 0. In quadratic form, this is
    # t² - Σᵢ xᵢ² ≥ 0 and t ≥ 0.

    # First, check the lower bound on t.

    _update_if_necessary(model)
    t_info = _info(model, f.variables[1])
    lb = _get_variable_lower_bound(model, t_info)
    if isnan(t_info.lower_bound_if_soc) && lb < 0.0
        # If `t_info.lower_bound_if_bounded` is active, this just makes
        # `t_info.lower_bound_if_soc` equal to it. If `lower_bound_if_bounded`
        # is set after, then it will call `_set_variable_lower_bound` and
        # update `lower_bound_if_soc` accordingly.
        t_info.lower_bound_if_soc = lb
        ret = GRBsetdblattrelement(model, "LB", Cint(t_info.column - 1), 0.0)
        _check_ret(model, ret)
    end
    t_info.num_soc_constraints += 1

    # Now add the quadratic constraint.

    I = Cint[column(model, v) - 1 for v in f.variables]
    V = fill(Cdouble(-1.0), length(f.variables))
    V[1] = 1.0
    ret = GRBaddqconstr(
        model,
        0,
        Cint[],
        Cdouble[],
        length(I),
        I,
        I,
        V,
        Cchar('>'),
        0.0,
        "",
    )
    _check_ret(model, ret)
    _require_update(model)
    model.last_constraint_index += 1
    model.quadratic_constraint_info[model.last_constraint_index] =
        _ConstraintInfo(length(model.quadratic_constraint_info) + 1, s)
    return MOI.ConstraintIndex{MOI.VectorOfVariables,MOI.SecondOrderCone}(
        model.last_constraint_index,
    )
end

function MOI.is_valid(
    model::Optimizer,
    c::MOI.ConstraintIndex{MOI.VectorOfVariables,MOI.SecondOrderCone},
)
    info = get(model.quadratic_constraint_info, c.value, nothing)
    return info !== nothing && typeof(info.set) == MOI.SecondOrderCone
end

function MOI.delete(
    model::Optimizer,
    c::MOI.ConstraintIndex{MOI.VectorOfVariables,MOI.SecondOrderCone},
)
    _update_if_necessary(model)
    f = MOI.get(model, MOI.ConstraintFunction(), c)
    info = _info(model, c)
    ret = GRBdelqconstrs(model, 1, [Cint(info.row - 1)])
    _check_ret(model, ret)
    _require_update(model)
    for (key, info_2) in model.quadratic_constraint_info
        if info_2.row > info.row
            info_2.row -= 1
        end
    end
    model.name_to_constraint_index = nothing
    delete!(model.quadratic_constraint_info, c.value)
    # Reset the lower bound on the `t` variable.
    t_info = _info(model, f.variables[1])
    t_info.num_soc_constraints -= 1
    if t_info.num_soc_constraints > 0
        # Don't do anything. There are still SOC associated with this variable.
        return
    elseif isnan(t_info.lower_bound_if_soc)
        # Don't do anything. It must have a >0 lower bound anyway.
        return
    end
    # There was a previous bound that we over-wrote, and it must have been
    # < 0 otherwise we wouldn't have needed to overwrite it.
    @assert t_info.lower_bound_if_soc < 0.0
    # Also, if there is a cached value in `t_info.lower_bound_if_bounded`
    # (i.e., `t_info.bound` is not `_NONE` nor `_LESS_THAN`), then it has
    # followed any changes `t_info.lower_bound_if_soc` has gone through
    # and has the same value. So when LB is set to the old value of
    # `lower_bound_if_soc` below, then `lower_bound_if_bounded` will stay
    # correct.
    @assert isnan(t_info.lower_bound_if_bounded) ||
            t_info.lower_bound_if_bounded == t_info.lower_bound_if_soc
    tmp_lower_bound = t_info.lower_bound_if_soc
    t_info.lower_bound_if_soc = NaN
    ret = GRBsetdblattrelement(
        model,
        "LB",
        Cint(t_info.column - 1),
        tmp_lower_bound,
    )
    _check_ret(model, ret)
    _require_update(model)
    return
end

function MOI.get(
    model::Optimizer,
    ::MOI.ConstraintSet,
    c::MOI.ConstraintIndex{MOI.VectorOfVariables,MOI.SecondOrderCone},
)
    return _info(model, c).set
end

function MOI.get(
    model::Optimizer,
    ::MOI.ConstraintFunction,
    c::MOI.ConstraintIndex{MOI.VectorOfVariables,MOI.SecondOrderCone},
)
    _update_if_necessary(model)
    a, b, I, J, V = _getqconstr(model, Cint(_info(model, c).row - 1))
    @assert length(a) == length(b) == 0  # Check for no linear terms.
    t = nothing
    x = MOI.VariableIndex[]
    for (i, j, coef) in zip(I, J, V)
        v = model.variable_info[CleverDicts.LinearIndex(i + 1)].index
        @assert i == j  # Check for no off-diagonals.
        if coef == 1.0
            @assert t === nothing  # There should only be one `t`.
            t = v
        else
            @assert coef == -1.0  # The coefficients _must_ be -1 for `x` terms.
            push!(x, v)
        end
    end
    @assert t !== nothing  # Check that we found a `t` variable.
    return MOI.VectorOfVariables([t; x])
end

function MOI.get(
    model::Optimizer,
    ::MOI.ConstraintPrimal,
    c::MOI.ConstraintIndex{MOI.VectorOfVariables,MOI.SecondOrderCone},
)
    f = MOI.get(model, MOI.ConstraintFunction(), c)
    return MOI.get(model, MOI.VariablePrimal(), f.variables)
end

function MOI.get(
    model::Optimizer,
    ::MOI.ConstraintName,
    c::MOI.ConstraintIndex{MOI.VectorOfVariables,MOI.SecondOrderCone},
)
    return _info(model, c).name
end

function MOI.set(
    model::Optimizer,
    ::MOI.ConstraintName,
    c::MOI.ConstraintIndex{MOI.VectorOfVariables,MOI.SecondOrderCone},
    name::String,
)
    info = _info(model, c)
    if !isempty(info.name) && model.name_to_constraint_index !== nothing
        delete!(model.name_to_constraint_index, info.name)
    end
    _update_if_necessary(model)
    ret = GRBsetstrattrelement(model, "QCName", Cint(info.row - 1), name)
    _check_ret(model, ret)
    _require_update(model)
    info.name = name
    if model.name_to_constraint_index === nothing || isempty(name)
        return
    end
    if haskey(model.name_to_constraint_index, name)
        model.name_to_constraint_index = nothing
    else
        model.name_to_constraint_index[c] = name
    end
    return
end
