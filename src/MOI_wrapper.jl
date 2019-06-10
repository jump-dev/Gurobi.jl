# Attributes still to implement:
#  - RawParameter
#  - BasisStatusCode
#  - ConstraintBasisStatus

import MathOptInterface

const MOI  = MathOptInterface

const SUPPORTED_OBJECTIVES = (
    MOI.SingleVariable,
    MOI.ScalarAffineFunction{Float64},
    MOI.ScalarQuadraticFunction{Float64}
)

const SUPPORTED_CONSTRAINTS = (
    (MOI.SingleVariable,                    MOI.EqualTo{Float64}),
    (MOI.SingleVariable,                    MOI.LessThan{Float64}),
    (MOI.SingleVariable,                    MOI.GreaterThan{Float64}),
    (MOI.SingleVariable,                    MOI.Interval{Float64}),

    (MOI.SingleVariable,                    MOI.ZeroOne),
    (MOI.SingleVariable,                    MOI.Integer),
    (MOI.SingleVariable,                    MOI.Semicontinuous{Float64}),
    (MOI.SingleVariable,                    MOI.Semiinteger{Float64}),

    (MOI.VectorOfVariables,                 MOI.SOS1{Float64}),
    (MOI.VectorOfVariables,                 MOI.SOS2{Float64}),

    (MOI.ScalarAffineFunction{Float64},     MOI.EqualTo{Float64}),
    (MOI.ScalarAffineFunction{Float64},     MOI.LessThan{Float64}),
    (MOI.ScalarAffineFunction{Float64},     MOI.GreaterThan{Float64}),
    # We choose _not_ to support ScalarAffineFunction-in-Interval because Gurobi
    # introduces some slack variables that makes it hard to keep track.
    # (MOI.ScalarAffineFunction{Float64},   MOI.Interval{Float64}),

    (MOI.ScalarQuadraticFunction{Float64},  MOI.EqualTo{Float64}),
    (MOI.ScalarQuadraticFunction{Float64},  MOI.LessThan{Float64}),
    (MOI.ScalarQuadraticFunction{Float64},  MOI.GreaterThan{Float64})
)

const SCALAR_SETS = Union{
    MOI.GreaterThan{Float64}, MOI.LessThan{Float64},
    MOI.EqualTo{Float64}, MOI.Interval{Float64}
}

@enum(VariableType, CONTINUOUS, BINARY, INTEGER, SEMIINTEGER, SEMICONTINUOUS)
@enum(BoundType, NONE, LESS_THAN, GREATER_THAN, LESS_AND_GREATER_THAN, INTERVAL, EQUAL_TO)
@enum(ObjectiveType, SINGLE_VARIABLE, SCALAR_AFFINE, SCALAR_QUADRATIC)

# Annoyingly, Gurobi ignores the name `""`, so we need a new default name.
const DEFAULT_NAME = " "

mutable struct VariableInfo
    column::Int
    bound::BoundType
    type::VariableType
    name::String
    # Storage for constraint names associated with variables because Gurobi
    # can only store names for variables and proper constraints.
    lessthan_name::String
    greaterthan_name::String
    equalto_name::String
    interval_name::String
    zeroone_name::String
    integer_name::String
    semiinteger_name::String
    semicontinuous_name::String
    VariableInfo(column::Int) = new(column, NONE, CONTINUOUS, "", "", "", "", "", "", "", "", "")
end

mutable struct ConstraintInfo
    row::Int
    set::MOI.AbstractSet
    # Storage for constraint names. Where possible, these are also stored in the
    # Gurobi model.
    name::String
    ConstraintInfo(row::Int, set) = new(row, set, "")
end

mutable struct Optimizer <: MOI.ModelLike
    inner::Model
    env::Union{Nothing, Env}
    silent::Bool
    params::Dict{String, Any}

    # The next field is used to cleverly manage calls to `update_model!`.
    # `needs_update` is used to record whether an update should be called before
    # accessing a model attribute (such as the value of a RHS term).
    needs_update::Bool

    objective_type::ObjectiveType
    is_feasibility::Bool

    last_variable_index::Int
    variable_info::Dict{MOI.VariableIndex, VariableInfo}
    columns::Vector{MOI.VariableIndex}

    last_constraint_index::Int
    affine_constraint_info::Dict{Int, ConstraintInfo}
    quadratic_constraint_info::Dict{Int, ConstraintInfo}
    sos_constraint_info::Dict{Int, ConstraintInfo}

    name_to_variable::Union{Nothing, Dict{String, MOI.VariableIndex}}
    name_to_constraint_index::Union{Nothing, Dict{String, MOI.ConstraintIndex}}

    has_unbounded_ray::Bool
    has_infeasibility_cert::Bool

    callback_variable_primal::Vector{Float64}

    """
        Optimizer(env = nothing; kwargs...)

    Create a new Optimizer object.

    You can share Gurobi `Env`s between models by passing an instance of `Env`
    as the first argument. By default, a new environment is created for every
    model.

    Note that we set the parameter `InfUnbdInfo` to `1` rather than the default
    of `0` so that we can query infeasibility certificates. Users are, however,
    free to overide this as follows `Gurobi.Optimizer(InfUndbInfo=0)`.
    """
    function Optimizer(env::Union{Nothing, Env} = nothing; kwargs...)
        model = new()
        model.env = env
        model.silent = false
        model.params = Dict{String, Any}()
        model.variable_info = Dict{MOI.VariableIndex, VariableInfo}()
        model.columns = MOI.VariableIndex[]
        model.affine_constraint_info = Dict{Int, ConstraintInfo}()
        model.quadratic_constraint_info = Dict{Int, Int}()
        model.sos_constraint_info = Dict{Int, ConstraintInfo}()
        model.last_variable_index = 0
        model.last_constraint_index = 0
        model.callback_variable_primal = Float64[]
        MOI.empty!(model)  # MOI.empty!(model) re-sets the `.inner` field.
        for (name, value) in kwargs
            model.params[string(name)] = value
            setparam!(model.inner, string(name), value)
        end
        if !haskey(model.params, "InfUnbdInfo")
            model.params["InfUnbdInfo"] = 1
            setparam!(model.inner, "InfUnbdInfo", 1)
        end
        return model
    end
end

function MOI.empty!(model::Optimizer)
    if model.env === nothing
        model.inner = Model(Env(), "", finalize_env = true)
    else
        model.inner = Model(model.env, "", finalize_env = false)
    end
    for (name, value) in model.params
        setparam!(model.inner, name, value)
    end
    if model.silent
        setparam!(model.inner, "OutputFlag", 0)
    end
    model.needs_update = false
    model.objective_type = SCALAR_AFFINE
    model.is_feasibility = true
    empty!(model.variable_info)
    model.name_to_variable = nothing
    empty!(model.columns)
    empty!(model.affine_constraint_info)
    model.name_to_constraint_index = nothing
    empty!(model.quadratic_constraint_info)
    empty!(model.sos_constraint_info)
    model.has_unbounded_ray = false
    model.has_infeasibility_cert = false
    empty!(model.callback_variable_primal)
    return
end

function MOI.is_empty(model::Optimizer)
    model.needs_update && return false
    model.objective_type != SCALAR_AFFINE && return false
    length(model.variable_info) != 0 && return false
    length(model.affine_constraint_info) != 0 && return false
    length(model.quadratic_constraint_info) != 0 && return false
    length(model.sos_constraint_info) != 0 && return false
    model.name_to_variable !== nothing && return false
    model.name_to_constraint_index !== nothing && return false
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
    _require_update(model::Optimizer)

Calls `update_model!`, but only if the `model.needs_update` flag is set.
"""
function _update_if_necessary(model::Optimizer)
    if model.needs_update
        update_model!(model.inner)
        model.needs_update = false
    end
    return
end

MOI.get(::Optimizer, ::MOI.SolverName) = "Gurobi"

function MOI.supports_constraint(::Optimizer, F::Type{<:MOI.AbstractFunction}, S::Type{<:MOI.AbstractSet})
    return (F, S) in SUPPORTED_CONSTRAINTS
end

function MOI.supports(::Optimizer, ::MOI.ObjectiveFunction{F}) where {F}
    return F in SUPPORTED_OBJECTIVES
end

MOI.supports(::Optimizer, ::MOI.VariableName, ::Type{MOI.VariableIndex}) = true
MOI.supports(::Optimizer, ::MOI.ConstraintName, ::Type{<:MOI.ConstraintIndex}) = true
# TODO: this next one looks like it might be a bug in MOI.Test.
MOI.supports(::Optimizer, ::MOI.ConstraintName, ::MOI.ConstraintIndex) = true
MOI.supports(::Optimizer, ::MOI.ObjectiveFunctionType) = true

MOI.supports(::Optimizer, ::MOI.Name) = true
MOI.supports(::Optimizer, ::MOI.Silent) = true
MOI.supports(::Optimizer, ::MOI.ConstraintSet, c) = true
MOI.supports(::Optimizer, ::MOI.ConstraintFunction, c) = true
MOI.supports(::Optimizer, ::MOI.ConstraintPrimal, c) = true
MOI.supports(::Optimizer, ::MOI.ConstraintDual, c) = true
MOI.supports(::Optimizer, ::MOI.ObjectiveSense) = true
MOI.supports(::Optimizer, ::MOI.ListOfConstraintIndices) = true
MOI.supports(::Optimizer, ::MOI.RawStatusString) = true

MOI.Utilities.supports_default_copy_to(::Optimizer, ::Bool) = true

function MOI.copy_to(dest::Optimizer, src::MOI.ModelLike; kwargs...)
    return MOI.Utilities.automatic_copy_to(dest, src; kwargs...)
end

function MOI.get(model::Optimizer, ::MOI.ListOfVariableAttributesSet)
    return MOI.AbstractVariableAttribute[MOI.VariableName()]
end

function MOI.get(model::Optimizer, ::MOI.ListOfModelAttributesSet)
    attributes = [
        MOI.ObjectiveSense(),
        MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}()
    ]
    if MOI.get(model, MOI.Name()) != ""
        push!(attributes, MOI.Name())
    end
    return attributes
end

function MOI.get(model::Optimizer, ::MOI.ListOfConstraintAttributesSet)
    return MOI.AbstractConstraintAttribute[MOI.ConstraintName()]
end

# These strings are taken directly from the following page of the online Gurobi
# documentation: https://www.gurobi.com/documentation/8.1/refman/optimization_status_codes.html#sec:StatusCodes
const RAW_STATUS_STRINGS = [
    "Model is loaded, but no solution information is available.",
    "Model was solved to optimality (subject to tolerances), and an optimal solution is available.",
    "Model was proven to be infeasible.",
    "Model was proven to be either infeasible or unbounded. To obtain a more definitive conclusion, set the DualReductions parameter to 0 and reoptimize.",
    "Model was proven to be unbounded. Important note: an unbounded status indicates the presence of an unbounded ray that allows the objective to improve without limit. It says nothing about whether the model has a feasible solution. If you require information on feasibility, you should set the objective to zero and reoptimize.",
    "Optimal objective for model was proven to be worse than the value specified in the Cutoff parameter. No solution information is available.",
    "Optimization terminated because the total number of simplex iterations performed exceeded the value specified in the IterationLimit parameter, or because the total number of barrier iterations exceeded the value specified in the BarIterLimit parameter.",
    "Optimization terminated because the total number of branch-and-cut nodes explored exceeded the value specified in the NodeLimit parameter.",
    "Optimization terminated because the time expended exceeded the value specified in the TimeLimit parameter.",
    "Optimization terminated because the number of solutions found reached the value specified in the SolutionLimit parameter.",
    "Optimization was terminated by the user.",
    "Optimization was terminated due to unrecoverable numerical difficulties.",
    "Unable to satisfy optimality tolerances; a sub-optimal solution is available.",
    "An asynchronous optimization call was made, but the associated optimization run is not yet complete.",
    "User specified an objective limit (a bound on either the best objective or the best bound), and that limit has been reached."
]

function MOI.get(model::Optimizer, ::MOI.RawStatusString)
    status_code = get_status_code(model.inner)
    return RAW_STATUS_STRINGS[status_code]
end

function indices_and_coefficients(model::Optimizer, f::MOI.ScalarAffineFunction{Float64})
    f_canon = MOI.Utilities.canonical(f)
    indices = Int[]
    coefficients = Float64[]
    for term in f_canon.terms
        push!(indices, model[term.variable_index].column)
        push!(coefficients, term.coefficient)
    end
    return indices, coefficients
end

sense_and_rhs(s::MOI.LessThan{Float64}) = (Cchar('<'), s.upper)
sense_and_rhs(s::MOI.GreaterThan{Float64}) = (Cchar('>'), s.lower)
sense_and_rhs(s::MOI.EqualTo{Float64}) = (Cchar('='), s.value)

function indices_and_coefficients(
    model::Optimizer, f::MOI.ScalarQuadraticFunction
)
    f_canon = MOI.Utilities.canonical(f)
    I, J, V = Int[], Int[], Float64[]
    for term in f_canon.quadratic_terms
        push!(I, model[term.variable_index_1].column)
        push!(J, model[term.variable_index_2].column)
        push!(V, term.coefficient)
        if I[end] == J[end]
            V[end] *= 0.5
        end
    end
    indices, coefficients = Int[], Float64[]
    for term in f_canon.affine_terms
        push!(indices, model[term.variable_index].column)
        push!(coefficients, term.coefficient)
    end
    return indices, coefficients, I, J, V
end

###
### Variables
###

# Short-cuts to return the VariableInfo associated with an index.
function Base.getindex(model::Optimizer, key::MOI.VariableIndex)
    if haskey(model.variable_info, key)
        return model.variable_info[key]
    end
    throw(MOI.InvalidIndex(key))
end

function MOI.add_variable(model::Optimizer)
    model.last_variable_index += 1
    index = MOI.VariableIndex(model.last_variable_index)
    model.variable_info[index] = VariableInfo(length(model.variable_info) + 1)
    add_cvar!(model.inner, 0.0)
    _require_update(model)
    push!(model.columns, index)
    return index
end

function MOI.add_variables(model::Optimizer, N::Int)
    Gurobi.add_cvars!(model.inner, zeros(N))
    _require_update(model)
    indices = Vector{MOI.VariableIndex}(undef, N)
    num_variables = length(model.variable_info)
    for i in 1:N
        model.last_variable_index += 1
        index = MOI.VariableIndex(model.last_variable_index)
        model.variable_info[index] = VariableInfo(num_variables + i)
        indices[i] = index
        push!(model.columns, index)
    end
    return indices
end

MOI.is_valid(model::Optimizer, v::MOI.VariableIndex) = haskey(model.variable_info, v)

function MOI.delete(model::Optimizer, v::MOI.VariableIndex)
    _update_if_necessary(model)
    info = model[v]
    del_vars!(model.inner, info.column)
    _require_update(model)
    for other_info in values(model.variable_info)
        if other_info.column > info.column
            other_info.column -= 1
        end
    end
    splice!(model.columns, info.column)
    model.name_to_variable = nothing
    delete!(model.variable_info, v)
    return
end

function MOI.get(model::Optimizer, ::Type{MOI.VariableIndex}, name::String)
    if model.name_to_variable === nothing
        rebuild_name_to_variable(model)
    end
    return get(model.name_to_variable, name, nothing)
end

function rebuild_name_to_variable(model::Optimizer)
    model.name_to_variable = Dict{String, MOI.VariableIndex}()
    for (index, info) in model.variable_info
        if info.name == ""
            continue
        end
        if haskey(model.name_to_variable, info.name)
            model.name_to_variable = nothing
            error("Duplicate variable name detected: $(info.name)")
        end
        model.name_to_variable[info.name] = index
    end
    return
end

function MOI.get(model::Optimizer, ::MOI.VariableName, v::MOI.VariableIndex)
    return model[v].name
end

function MOI.set(model::Optimizer, ::MOI.VariableName, v::MOI.VariableIndex, name::String)
    info = model[v]
    info.name = name
    _update_if_necessary(model)
    set_strattrelement!(model.inner, "VarName", info.column, name == "" ? " " : name)
    _require_update(model)
    if model.name_to_variable !== nothing && !haskey(model.name_to_variable, name)
        model.name_to_variable[name] = v
    end
    return
end

###
### Objectives
###

function MOI.set(
    model::Optimizer, ::MOI.ObjectiveSense, sense::MOI.OptimizationSense
)
    if sense == MOI.MIN_SENSE
        set_sense!(model.inner, :minimize)
        model.is_feasibility = false
    elseif sense == MOI.MAX_SENSE
        set_sense!(model.inner, :maximize)
        model.is_feasibility = false
    elseif sense == MOI.FEASIBILITY_SENSE
        set_sense!(model.inner, :minimize)
        model.is_feasibility = true
    else
        error("Invalid objective sense: $(sense)")
    end
    _require_update(model)
    return
end

function MOI.get(model::Optimizer, ::MOI.ObjectiveSense)
    _update_if_necessary(model)
    sense = model_sense(model.inner)
    if model.is_feasibility
        return MOI.FEASIBILITY_SENSE
    elseif sense == :maximize
        return MOI.MAX_SENSE
    elseif sense == :minimize
        return MOI.MIN_SENSE
    end
    error("Invalid objective sense: $(sense)")
end

function MOI.set(
    model::Optimizer, ::MOI.ObjectiveFunction{F}, f::F
) where {F <: MOI.SingleVariable}
    MOI.set(
        model, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(),
        convert(MOI.ScalarAffineFunction{Float64}, f)
    )
    model.objective_type = SINGLE_VARIABLE
    return
end

function MOI.get(model::Optimizer, ::MOI.ObjectiveFunction{MOI.SingleVariable})
    obj = MOI.get(model, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}())
    return convert(MOI.SingleVariable, obj)
end

function MOI.set(
    model::Optimizer, ::MOI.ObjectiveFunction{F}, f::F
) where {F <: MOI.ScalarAffineFunction{Float64}}
    if model.objective_type == SCALAR_QUADRATIC
        # We need to zero out the existing quadratic objective.
        delq!(model.inner)
    end
    num_vars = length(model.variable_info)
    obj = zeros(Float64, num_vars)
    for term in f.terms
        column = model[term.variable_index].column
        obj[column] += term.coefficient
    end
    # This update is needed because we might have added some variables.
    _update_if_necessary(model)
    set_dblattrarray!(model.inner, "Obj", 1, num_vars, obj)
    set_dblattr!(model.inner, "ObjCon", f.constant)
    _require_update(model)
    model.objective_type = SCALAR_AFFINE
end

function MOI.get(
    model::Optimizer, ::MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}
)
    if model.objective_type == SCALAR_QUADRATIC
        error("Unable to get objective function. Currently: $(model.objective_type).")
    end
    _update_if_necessary(model)
    dest = zeros(length(model.variable_info))
    get_dblattrarray!(dest, model.inner, "Obj", 1)
    terms = MOI.ScalarAffineTerm{Float64}[]
    for (index, info) in model.variable_info
        coefficient = dest[info.column]
        iszero(coefficient) && continue
        push!(terms, MOI.ScalarAffineTerm(coefficient, index))
    end
    constant = get_dblattr(model.inner, "ObjCon")
    return MOI.ScalarAffineFunction(terms, constant)
end

function MOI.set(
    model::Optimizer, ::MOI.ObjectiveFunction{F}, f::F
) where {F <: MOI.ScalarQuadraticFunction{Float64}}
    affine_indices, affine_coefficients, I, J, V = indices_and_coefficients(model, f)
    _update_if_necessary(model)
    # We need to zero out any existing linear objective.
    obj = zeros(length(model.variable_info))
    for (i, c) in zip(affine_indices, affine_coefficients)
        obj[i] = c
    end
    set_dblattrarray!(model.inner, "Obj", 1, length(obj), obj)
    set_dblattr!(model.inner, "ObjCon", f.constant)
    # We need to zero out the existing quadratic objective.
    delq!(model.inner)
    add_qpterms!(model.inner, I, J, V)
    _require_update(model)
    model.objective_type = SCALAR_QUADRATIC
    return
end

function MOI.get(
    model::Optimizer,
    ::MOI.ObjectiveFunction{MOI.ScalarQuadraticFunction{Float64}}
)
    _update_if_necessary(model)
    dest = zeros(length(model.variable_info))
    get_dblattrarray!(dest, model.inner, "Obj", 1)
    terms = MOI.ScalarAffineTerm{Float64}[]
    for (index, info) in model.variable_info
        coefficient = dest[info.column]
        iszero(coefficient) && continue
        push!(terms, MOI.ScalarAffineTerm(coefficient, index))
    end
    constant = get_dblattr(model.inner, "ObjCon")
    q_terms = MOI.ScalarQuadraticTerm{Float64}[]
    # getq returns a list of terms. MOI requires 0.5 x' Q x. So, to get from
    #   Gurobi -> MOI => multiply diagonals by 2.0
    #   MOI -> Gurobi => multiply diagonals by 0.5
    # Example: 2x^2 + x*y + y^2
    #   Gurobi returns: (I, J, V) = ([0, 0, 1], [0, 1, 1], [2, 1, 1])
    #   MOI needs:
    #     [SQT(4.0, x, x), SQT(1.0, x, y), SQT(2.0, y, y)]
    #   |x y| * |a b| * |x| = |ax+by bx+cy| * |x| = 0.5ax^2 + bxy + 0.5cy^2
    #           |b c|   |y|                   |y|
    I, J, V = getq(model.inner)
    for (i, j, v) in zip(I, J, V)
        iszero(v) && continue
        new_v = i == j ? 2v : v
        push!(
            q_terms,
            MOI.ScalarQuadraticTerm(
                new_v, model.columns[i + 1], model.columns[j + 1]
            )
        )
    end
    return MOI.ScalarQuadraticFunction(terms, q_terms, constant)
end

function MOI.modify(
    model::Optimizer,
    ::MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}},
    chg::MOI.ScalarConstantChange{Float64}
)
    set_dblattr!(model.inner, "ObjCon", chg.new_constant)
    _require_update(model)
    return
end

##
##  SingleVariable-in-Set constraints.
##

function Base.getindex(model::Optimizer, c::MOI.ConstraintIndex{MOI.SingleVariable, <:Any})
    var_index = MOI.VariableIndex(c.value)
    if haskey(model.variable_info, var_index)
        return model[var_index]
    end
    return throw(MOI.InvalidIndex(c))
end

function throw_if_invalid(model, c::MOI.ConstraintIndex{MOI.SingleVariable, <:Any})
    if !MOI.is_valid(model, c)
        throw(MOI.InvalidIndex(c))
    end
end

function MOI.is_valid(
    model::Optimizer,
    c::MOI.ConstraintIndex{MOI.SingleVariable, MOI.LessThan{Float64}}
)
    if haskey(model.variable_info, MOI.VariableIndex(c.value))
        info = model[c]
        return info.bound == LESS_THAN || info.bound == LESS_AND_GREATER_THAN
    end
    return false
end

function MOI.is_valid(
    model::Optimizer,
    c::MOI.ConstraintIndex{MOI.SingleVariable, MOI.GreaterThan{Float64}}
)
    if haskey(model.variable_info, MOI.VariableIndex(c.value))
        info = model[c]
        return info.bound == GREATER_THAN || info.bound == LESS_AND_GREATER_THAN
    end
    return false
end

function MOI.is_valid(
    model::Optimizer,
    c::MOI.ConstraintIndex{MOI.SingleVariable, MOI.Interval{Float64}}
)
    return haskey(model.variable_info, MOI.VariableIndex(c.value)) &&
        model[c].bound == INTERVAL
end

function MOI.is_valid(
    model::Optimizer,
    c::MOI.ConstraintIndex{MOI.SingleVariable, MOI.EqualTo{Float64}}
)
    return haskey(model.variable_info, MOI.VariableIndex(c.value)) &&
        model[c].bound == EQUAL_TO
end

function MOI.is_valid(
    model::Optimizer,
    c::MOI.ConstraintIndex{MOI.SingleVariable, MOI.ZeroOne}
)
    return haskey(model.variable_info, MOI.VariableIndex(c.value)) &&
        model[c].type == BINARY
end

function MOI.is_valid(
    model::Optimizer,
    c::MOI.ConstraintIndex{MOI.SingleVariable, MOI.Integer}
)
    return haskey(model.variable_info, MOI.VariableIndex(c.value)) &&
        model[c].type == INTEGER
end

function MOI.is_valid(
    model::Optimizer,
    c::MOI.ConstraintIndex{MOI.SingleVariable, MOI.Semicontinuous{Float64}}
)
    return haskey(model.variable_info, MOI.VariableIndex(c.value)) &&
        model[c].type == SEMICONTINUOUS
end

function MOI.is_valid(
    model::Optimizer,
    c::MOI.ConstraintIndex{MOI.SingleVariable, MOI.Semiinteger{Float64}}
)
    return haskey(model.variable_info, MOI.VariableIndex(c.value)) &&
        model[c].type == SEMIINTEGER
end

function MOI.get(
    model::Optimizer, ::MOI.ConstraintFunction,
    c::MOI.ConstraintIndex{MOI.SingleVariable, <:Any}
)
    throw_if_invalid(model, c)
    return MOI.SingleVariable(MOI.VariableIndex(c.value))
end

function MOI.set(
    model::Optimizer, ::MOI.ConstraintFunction,
    c::MOI.ConstraintIndex{MOI.SingleVariable, <:Any}, ::MOI.SingleVariable
)
    return throw(MOI.SettingSingleVariableFunctionNotAllowed())
end

bounds(s::MOI.GreaterThan{Float64}) = (s.lower, nothing)
bounds(s::MOI.LessThan{Float64}) = (nothing, s.upper)
bounds(s::MOI.EqualTo{Float64}) = (s.value, s.value)
bounds(s::MOI.Interval{Float64}) = (s.lower, s.upper)

function throw_if_existing_lower(bound, var_type, new_set, variable)
    existing_set = if bound == LESS_AND_GREATER_THAN || bound == GREATER_THAN
        MOI.GreaterThan{Float64}
    elseif bound == INTERVAL
        MOI.Interval{Float64}
    elseif bound == EQUAL_TO
        MOI.EqualTo{Float64}
    elseif var_type == SEMIINTEGER
        MOI.Semiinteger{Float64}
    elseif var_type == SEMICONTINUOUS
        MOI.Semicontinuous{Float64}
    else
        nothing  # Also covers `NONE` and `LESS_THAN`.
    end
    if existing_set !== nothing
        throw(MOI.LowerBoundAlreadySet{existing_set, new_set}(variable))
    end
end

function throw_if_existing_upper(bound, var_type, new_set, variable)
    existing_set = if bound == LESS_AND_GREATER_THAN || bound == LESS_THAN
        MOI.LessThan{Float64}
    elseif bound == INTERVAL
        MOI.Interval{Float64}
    elseif bound == EQUAL_TO
        MOI.EqualTo{Float64}
    elseif var_type == SEMIINTEGER
        MOI.Semiinteger{Float64}
    elseif var_type == SEMICONTINUOUS
        MOI.Semicontinuous{Float64}
    else
        nothing  # Also covers `NONE` and `GREATER_THAN`.
    end
    if existing_set !== nothing
        throw(MOI.UpperBoundAlreadySet{existing_set, new_set}(variable))
    end
end

function MOI.add_constraint(
    model::Optimizer, f::MOI.SingleVariable, s::S
) where {S <: SCALAR_SETS}
    info = model[f.variable]
    if typeof(s) == MOI.LessThan{Float64}
        throw_if_existing_upper(info.bound, info.type, S, f.variable)
        info.bound = info.bound == GREATER_THAN ? LESS_AND_GREATER_THAN : LESS_THAN
    elseif typeof(s) == MOI.GreaterThan{Float64}
        throw_if_existing_lower(info.bound, info.type, S, f.variable)
        info.bound = info.bound == LESS_THAN ? LESS_AND_GREATER_THAN : GREATER_THAN
    elseif typeof(s) == MOI.EqualTo{Float64}
        throw_if_existing_lower(info.bound, info.type, S, f.variable)
        throw_if_existing_upper(info.bound, info.type, S, f.variable)
        info.bound = EQUAL_TO
    elseif typeof(s) == MOI.Interval{Float64}
        throw_if_existing_lower(info.bound, info.type, S, f.variable)
        throw_if_existing_upper(info.bound, info.type, S, f.variable)
        info.bound = INTERVAL
    end
    index = MOI.ConstraintIndex{MOI.SingleVariable, typeof(s)}(f.variable.value)
    MOI.set(model, MOI.ConstraintSet(), index, s)
    return index
end

function MOI.delete(
    model::Optimizer, c::MOI.ConstraintIndex{MOI.SingleVariable, MOI.LessThan{Float64}}
)
    throw_if_invalid(model, c)
    _update_if_necessary(model)
    info = model[c]
    Gurobi.set_dblattrelement!(model.inner, "UB", info.column, Inf)
    _require_update(model)
    if info.bound == LESS_AND_GREATER_THAN
        info.bound = GREATER_THAN
    else
        info.bound = NONE
    end
    return
end

function MOI.delete(
    model::Optimizer, c::MOI.ConstraintIndex{MOI.SingleVariable, MOI.GreaterThan{Float64}}
)
    throw_if_invalid(model, c)
    info = model[c]
    _update_if_necessary(model)
    set_dblattrelement!(model.inner, "LB", info.column, -Inf)
    _require_update(model)
    if info.bound == LESS_AND_GREATER_THAN
        info.bound = LESS_THAN
    else
        info.bound = NONE
    end
    return
end

function MOI.delete(
    model::Optimizer, c::MOI.ConstraintIndex{MOI.SingleVariable, MOI.Interval{Float64}}
)
    throw_if_invalid(model, c)
    info = model[c]
    _update_if_necessary(model)
    Gurobi.set_dblattrelement!(model.inner, "LB", info.column, -Inf)
    Gurobi.set_dblattrelement!(model.inner, "UB", info.column, Inf)
    _require_update(model)
    info.bound = NONE
    return
end

function MOI.delete(
    model::Optimizer, c::MOI.ConstraintIndex{MOI.SingleVariable, MOI.EqualTo{Float64}}
)
    throw_if_invalid(model, c)
    info = model[c]
    _update_if_necessary(model)
    Gurobi.set_dblattrelement!(model.inner, "LB", info.column, -Inf)
    Gurobi.set_dblattrelement!(model.inner, "UB", info.column, Inf)
    _require_update(model)
    info.bound = NONE
    return
end

function MOI.get(
    model::Optimizer, ::MOI.ConstraintSet,
    c::MOI.ConstraintIndex{MOI.SingleVariable, MOI.GreaterThan{Float64}}
)
    throw_if_invalid(model, c)
    _update_if_necessary(model)
    lower = Gurobi.get_dblattrelement(model.inner, "LB", model[c].column)
    return MOI.GreaterThan(lower)
end

function MOI.get(
    model::Optimizer, ::MOI.ConstraintSet,
    c::MOI.ConstraintIndex{MOI.SingleVariable, MOI.LessThan{Float64}}
)
    throw_if_invalid(model, c)
    _update_if_necessary(model)
    upper = Gurobi.get_dblattrelement(model.inner, "UB", model[c].column)
    return MOI.LessThan(upper)
end

function MOI.get(
    model::Optimizer, ::MOI.ConstraintSet,
    c::MOI.ConstraintIndex{MOI.SingleVariable, MOI.EqualTo{Float64}}
)
    throw_if_invalid(model, c)
    _update_if_necessary(model)
    lower = Gurobi.get_dblattrelement(model.inner, "LB", model[c].column)
    return MOI.EqualTo(lower)
end

function MOI.get(
    model::Optimizer, ::MOI.ConstraintSet,
    c::MOI.ConstraintIndex{MOI.SingleVariable, MOI.Interval{Float64}}
)
    throw_if_invalid(model, c)
    _update_if_necessary(model)
    lower = Gurobi.get_dblattrelement(model.inner, "LB", model[c].column)
    upper = Gurobi.get_dblattrelement(model.inner, "UB", model[c].column)
    return MOI.Interval(lower, upper)
end

function MOI.set(
    model::Optimizer, ::MOI.ConstraintSet,
    c::MOI.ConstraintIndex{MOI.SingleVariable, S}, s::S
) where {S<:SCALAR_SETS}
    throw_if_invalid(model, c)
    lower, upper = bounds(s)
    info = model[c]
    _update_if_necessary(model)
    if lower !== nothing
        Gurobi.set_dblattrelement!(model.inner, "LB", info.column, lower)
    end
    if upper !== nothing
        Gurobi.set_dblattrelement!(model.inner, "UB", info.column, upper)
    end
    _require_update(model)
    return
end

function MOI.add_constraint(model::Optimizer, f::MOI.SingleVariable, ::MOI.ZeroOne)
    info = model[f.variable]
    _update_if_necessary(model)
    Gurobi.set_charattrelement!(model.inner, "VType", info.column, Char('B'))
    _require_update(model)
    info.type = BINARY
    return MOI.ConstraintIndex{MOI.SingleVariable, MOI.ZeroOne}(f.variable.value)
end

function MOI.delete(model::Optimizer, c::MOI.ConstraintIndex{MOI.SingleVariable, MOI.ZeroOne})
    throw_if_invalid(model, c)
    info = model[c]
    _update_if_necessary(model)
    Gurobi.set_charattrelement!(model.inner, "VType", info.column, Char('C'))
    _require_update(model)
    info.type = CONTINUOUS
    return
end

function MOI.get(
    model::Optimizer, ::MOI.ConstraintSet,
    c::MOI.ConstraintIndex{MOI.SingleVariable, MOI.ZeroOne}
)
    throw_if_invalid(model, c)
    return MOI.ZeroOne()
end

function MOI.add_constraint(
    model::Optimizer, f::MOI.SingleVariable, ::MOI.Integer
)
    info = model[f.variable]
    _update_if_necessary(model)
    Gurobi.set_charattrelement!(model.inner, "VType", info.column, Char('I'))
    _require_update(model)
    info.type = INTEGER
    return MOI.ConstraintIndex{MOI.SingleVariable, MOI.Integer}(f.variable.value)
end

function MOI.delete(
    model::Optimizer, c::MOI.ConstraintIndex{MOI.SingleVariable, MOI.Integer}
)
    throw_if_invalid(model, c)
    info = model[c]
    _update_if_necessary(model)
    Gurobi.set_charattrelement!(model.inner, "VType", info.column, Char('C'))
    _require_update(model)
    info.type = CONTINUOUS
    return
end

function MOI.get(
    model::Optimizer, ::MOI.ConstraintSet,
    c::MOI.ConstraintIndex{MOI.SingleVariable, MOI.Integer}
)
    throw_if_invalid(model, c)
    return MOI.Integer()
end

function MOI.add_constraint(
    model::Optimizer, f::MOI.SingleVariable, s::MOI.Semicontinuous{Float64}
)
    info = model[f.variable]
    throw_if_existing_lower(info.bound, info.type, typeof(s), f.variable)
    throw_if_existing_upper(info.bound, info.type, typeof(s), f.variable)
    _update_if_necessary(model)
    Gurobi.set_charattrelement!(model.inner, "VType", info.column, Char('S'))
    Gurobi.set_dblattrelement!(model.inner, "LB", info.column, s.lower)
    Gurobi.set_dblattrelement!(model.inner, "UB", info.column, s.upper)
    _require_update(model)
    info.type = SEMICONTINUOUS
    return MOI.ConstraintIndex{MOI.SingleVariable, MOI.Semicontinuous{Float64}}(f.variable.value)
end

function MOI.delete(
    model::Optimizer,
    c::MOI.ConstraintIndex{MOI.SingleVariable, MOI.Semicontinuous{Float64}}
)
    throw_if_invalid(model, c)
    info = model[c]
    _update_if_necessary(model)
    Gurobi.set_charattrelement!(model.inner, "VType", info.column, Char('C'))
    Gurobi.set_dblattrelement!(model.inner, "LB", info.column, -Inf)
    Gurobi.set_dblattrelement!(model.inner, "UB", info.column, Inf)
    _require_update(model)
    info.type = CONTINUOUS
    return
end

function MOI.get(
    model::Optimizer, ::MOI.ConstraintSet,
    c::MOI.ConstraintIndex{MOI.SingleVariable, MOI.Semicontinuous{Float64}}
)
    throw_if_invalid(model, c)
    info = model[c]
    _update_if_necessary(model)
    lower = Gurobi.get_dblattrelement(model.inner, "LB", info.column)
    upper = Gurobi.get_dblattrelement(model.inner, "UB", info.column)
    return MOI.Semicontinuous(lower, upper)
end

function MOI.add_constraint(
    model::Optimizer, f::MOI.SingleVariable, s::MOI.Semiinteger{Float64}
)
    info = model[f.variable]
    throw_if_existing_lower(info.bound, info.type, typeof(s), f.variable)
    throw_if_existing_upper(info.bound, info.type, typeof(s), f.variable)
    _update_if_necessary(model)
    Gurobi.set_charattrelement!(model.inner, "VType", info.column, Char('N'))
    Gurobi.set_dblattrelement!(model.inner, "LB", info.column, s.lower)
    Gurobi.set_dblattrelement!(model.inner, "UB", info.column, s.upper)
    _require_update(model)
    info.type = SEMIINTEGER
    return MOI.ConstraintIndex{MOI.SingleVariable, MOI.Semiinteger{Float64}}(f.variable.value)
end

function MOI.delete(
    model::Optimizer,
    c::MOI.ConstraintIndex{MOI.SingleVariable, MOI.Semiinteger{Float64}}
)
    throw_if_invalid(model, c)
    info = model[c]
    _update_if_necessary(model)
    Gurobi.set_charattrelement!(model.inner, "VType", info.column, Char('C'))
    Gurobi.set_dblattrelement!(model.inner, "LB", info.column, -Inf)
    Gurobi.set_dblattrelement!(model.inner, "UB", info.column, Inf)
    _require_update(model)
    info.type = CONTINUOUS
    return
end

function MOI.get(
    model::Optimizer, ::MOI.ConstraintSet,
    c::MOI.ConstraintIndex{MOI.SingleVariable, MOI.Semiinteger{Float64}}
)
    throw_if_invalid(model, c)
    info = model[c]
    _update_if_necessary(model)
    lower = Gurobi.get_dblattrelement(model.inner, "LB", info.column)
    upper = Gurobi.get_dblattrelement(model.inner, "UB", info.column)
    return MOI.Semiinteger(lower, upper)
end

constraint_name(info::VariableInfo, ::Type{<:MOI.LessThan}) = info.lessthan_name
constraint_name(info::VariableInfo, ::Type{<:MOI.GreaterThan}) = info.greaterthan_name
constraint_name(info::VariableInfo, ::Type{<:MOI.EqualTo}) = info.equalto_name
constraint_name(info::VariableInfo, ::Type{<:MOI.Interval}) = info.interval_name
constraint_name(info::VariableInfo, ::Type{<:MOI.ZeroOne}) = info.zeroone_name
constraint_name(info::VariableInfo, ::Type{<:MOI.Integer}) = info.integer_name
constraint_name(info::VariableInfo, ::Type{<:MOI.Semiinteger}) = info.semiinteger_name
constraint_name(info::VariableInfo, ::Type{<:MOI.Semicontinuous}) = info.semicontinuous_name

function MOI.get(
    model::Optimizer, ::MOI.ConstraintName,
    c::MOI.ConstraintIndex{MOI.SingleVariable, S}
) where {S}
    throw_if_invalid(model, c)
    return constraint_name(model[c], S)
end

function MOI.set(
    model::Optimizer, ::MOI.ConstraintName,
    c::MOI.ConstraintIndex{MOI.SingleVariable, S}, name::String
) where {S}
    throw_if_invalid(model, c)
    if S <: MOI.LessThan
        model[c].lessthan_name = name
    elseif S <: MOI.GreaterThan
        model[c].greaterthan_name = name
    elseif S <: MOI.Interval
        model[c].interval_name = name
    elseif S <: MOI.EqualTo
        model[c].equalto_name = name
    elseif S <: MOI.ZeroOne
        model[c].zeroone_name = name
    elseif S <: MOI.Integer
        model[c].integer_name = name
    elseif S <: MOI.Semicontinuous
        model[c].semicontinuous_name = name
    elseif S <: MOI.Semiinteger
        model[c].semiinteger_name = name
    end
    return
end


###
### ScalarAffineFunction-in-Set
###

function Base.getindex(
    model::Optimizer,
    key::MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64}, <:Any}
)
    if haskey(model.affine_constraint_info, key.value)
        return model.affine_constraint_info[key.value]
    end
    throw(MOI.InvalidIndex(key))
end

function MOI.is_valid(
    model::Optimizer,
    c::MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64}, S}
) where {S}
    info = get(model.affine_constraint_info, c.value, nothing)
    if info === nothing
        return false
    else
        return typeof(info.set) == S
    end
end

function MOI.add_constraint(
    model::Optimizer, f::MOI.ScalarAffineFunction{Float64},
    s::Union{MOI.GreaterThan{Float64}, MOI.LessThan{Float64}, MOI.EqualTo{Float64}}
)
    if !iszero(f.constant)
        throw(MOI.ScalarFunctionConstantNotZero{Float64, typeof(f), typeof(s)}(f.constant))
    end
    model.last_constraint_index += 1
    model.affine_constraint_info[model.last_constraint_index] =
        ConstraintInfo(length(model.affine_constraint_info) + 1, s)

    indices, coefficients = indices_and_coefficients(model, f)
    sense, rhs = sense_and_rhs(s)
    add_constr!(model.inner, indices, coefficients, sense, rhs)
    _require_update(model)
    return MOI.ConstraintIndex{typeof(f), typeof(s)}(model.last_constraint_index)
end

function MOI.delete(
    model::Optimizer,
    c::MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64}, <:Any}
)
    _update_if_necessary(model)
    del_constrs!(model.inner, model[c].row)
    _require_update(model)
    row = model[c].row
    for (key, info) in model.affine_constraint_info
        if info.row > row
            info.row -= 1
        end
    end
    model.name_to_constraint_index = nothing
    delete!(model.affine_constraint_info, c.value)
    return
end

function MOI.get(
    model::Optimizer, ::MOI.ConstraintSet,
    c::MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64}, S}
) where {S}
    _update_if_necessary(model)
    rhs = get_dblattrelement(model.inner, "RHS", model[c].row)
    return S(rhs)
end

function MOI.set(
    model::Optimizer, ::MOI.ConstraintSet,
    c::MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64}, S}, s::S
) where {S}
    _update_if_necessary(model)
    set_dblattrelement!(model.inner, "RHS", model[c].row, MOI.constant(s))
    _require_update(model)
    return
end

function MOI.get(
    model::Optimizer, ::MOI.ConstraintFunction,
    c::MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64}, S}
) where {S}
    _update_if_necessary(model)
    sparse_a = SparseArrays.sparse(get_constrs(model.inner, model[c].row, 1)')
    terms = MOI.ScalarAffineTerm{Float64}[]
    for (col, val) in zip(sparse_a.rowval, sparse_a.nzval)
        iszero(val) && continue
        push!(terms, MOI.ScalarAffineTerm(val, model.columns[col]))
    end
    return MOI.ScalarAffineFunction(terms, 0.0)
end

function MOI.get(
    model::Optimizer, ::MOI.ConstraintName,
    c::MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64}, <:Any}
)
    return model[c].name
end

function MOI.set(
    model::Optimizer, ::MOI.ConstraintName,
    c::MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64}, <:Any}, name::String
)
    model[c].name = name
    _update_if_necessary(model)
    set_strattrelement!(model.inner, "ConstrName", model[c].row, name == "" ? " " : name)
    _require_update(model)
    if model.name_to_constraint_index !== nothing && !haskey(model.name_to_constraint_index, name)
        model.name_to_constraint_index[c] = name
    else
        model.name_to_constraint_index = nothing
    end
    return
end

function MOI.get(model::Optimizer, ::Type{MOI.ConstraintIndex}, name::String)
    if model.name_to_constraint_index === nothing
        rebuild_name_to_constraint_index(model)
    end
    index = get(model.name_to_constraint_index, name, nothing)
    if index === nothing
        for S in (MOI.LessThan{Float64}, MOI.GreaterThan{Float64},
            MOI.EqualTo{Float64}, MOI.Interval{Float64}, MOI.ZeroOne,
            MOI.Integer, MOI.Semicontinuous{Float64}, MOI.Semiinteger{Float64}
        )
            index_2 = MOI.get(
                model, MOI.ConstraintIndex{MOI.SingleVariable, S}, name
            )
            if index_2 === nothing
                continue
            else
                if index === nothing
                    index = index_2
                else
                    error("Duplicate name detected: ", name)
                end
            end
        end
    end
    return index
end

function MOI.get(
    model::Optimizer, C::Type{MOI.ConstraintIndex{MOI.SingleVariable, S}}, name::String
) where {S}
    index = nothing
    for (key, info) in model.variable_info
        if constraint_name(info, S) == name
            if index === nothing
                index = key
            else
                error("Duplicate name detected: ", name)
            end
        end
    end
    if index === nothing
        return nothing
    end
    return MOI.ConstraintIndex{MOI.SingleVariable, S}(index.value)
end

function MOI.get(
    model::Optimizer, C::Type{MOI.ConstraintIndex{F, S}}, name::String
) where {F, S}
    index = MOI.get(model, MOI.ConstraintIndex, name)
    if typeof(index) == C
        return index::MOI.ConstraintIndex{F, S}
    end
    return nothing
end

function rebuild_name_to_constraint_index(model::Optimizer)
    model.name_to_constraint_index = Dict{String, Int}()
    _update_if_necessary(model)
    _rebuild_name_to_constraint_index_util(
        model, model.affine_constraint_info, MOI.ScalarAffineFunction{Float64}
    )
    _rebuild_name_to_constraint_index_util(
        model, model.quadratic_constraint_info, MOI.ScalarQuadraticFunction{Float64}
    )
    _rebuild_name_to_constraint_index_util(
        model, model.sos_constraint_info, MOI.VectorOfVariables
    )
    return
end

function _rebuild_name_to_constraint_index_util(model::Optimizer, dict, F)
    for (index, info) in dict
        info.name == "" && continue
        if haskey(model.name_to_constraint_index, info.name)
            model.name_to_constraint_index = nothing
            error("Duplicate variable name detected: $(info.name)")
        end
        model.name_to_constraint_index[info.name] =
            MOI.ConstraintIndex{F, typeof(info.set)}(index)
    end
    return
end

###
### ScalarQuadraticFunction-in-SCALAR_SET
###

function Base.getindex(
    model::Optimizer,
    c::MOI.ConstraintIndex{MOI.ScalarQuadraticFunction{Float64}, S}
) where {S}
    if haskey(model.quadratic_constraint_info, c.value)
        return model.quadratic_constraint_info[c.value]
    end
    throw(MOI.InvalidIndex(c))
end

function MOI.add_constraint(
    model::Optimizer,
    f::MOI.ScalarQuadraticFunction{Float64}, s::SCALAR_SETS
)
    if !iszero(f.constant)
        throw(MOI.ScalarFunctionConstantNotZero{Float64, typeof(f), typeof(s)}(f.constant))
    end
    indices, coefficients, I, J, V = indices_and_coefficients(model, f)
    sense, rhs = sense_and_rhs(s)
    _update_if_necessary(model)
    add_qconstr!(model.inner, indices, coefficients, I, J, V, sense, rhs)
    _update_if_necessary(model)
    _require_update(model)
    model.last_constraint_index += 1
    model.quadratic_constraint_info[model.last_constraint_index] = ConstraintInfo(
        length(model.quadratic_constraint_info) + 1, s
    )
    return MOI.ConstraintIndex{MOI.ScalarQuadraticFunction{Float64}, typeof(s)}(model.last_constraint_index)
end

function MOI.is_valid(
    model::Optimizer,
    c::MOI.ConstraintIndex{MOI.ScalarQuadraticFunction{Float64}, S}
) where {S}
    info = get(model.quadratic_constraint_info, c.value, nothing)
    return info !== nothing && typeof(info.set) == S
end

function MOI.delete(
    model::Optimizer,
    c::MOI.ConstraintIndex{MOI.ScalarQuadraticFunction{Float64}, S}
) where {S}
    _update_if_necessary(model)
    delqconstrs!(model.inner, [model[c].row])
    _require_update(model)
    row = model[c].row
    for (key, info) in model.quadratic_constraint_info
        if info.row > row
            info.row -= 1
        end
    end
    model.name_to_constraint_index = nothing
    delete!(model.quadratic_constraint_info, c.value)
    return
end

function MOI.get(
    model::Optimizer, ::MOI.ConstraintSet,
    c::MOI.ConstraintIndex{MOI.ScalarQuadraticFunction{Float64}, S}
) where {S}
    _update_if_necessary(model)
    rhs = get_dblattrelement(model.inner, "QCRHS", model[c].row)
    return S(rhs)
end

function MOI.get(
    model::Optimizer, ::MOI.ConstraintFunction,
    c::MOI.ConstraintIndex{MOI.ScalarQuadraticFunction{Float64}, S}
) where {S}
    _update_if_necessary(model)
    affine_cols, affine_coefficients, I, J, V = getqconstr(model.inner, model[c].row)
    affine_terms = MOI.ScalarAffineTerm{Float64}[]
    for (col, coef) in zip(affine_cols, affine_coefficients)
        iszero(coef) && continue
        push!(
            affine_terms,
            MOI.ScalarAffineTerm(coef, model.columns[col + 1])
        )
    end
    quadratic_terms = MOI.ScalarQuadraticTerm{Float64}[]
    for (i, j, coef) in zip(I, J, V)
        new_coef = i == j ? 2coef : coef
        push!(
            quadratic_terms,
            MOI.ScalarQuadraticTerm(new_coef, model.columns[i + 1], model.columns[j + 1])
        )
    end
    constant = get_dblattr(model.inner, "ObjCon")
    return MOI.ScalarQuadraticFunction(affine_terms, quadratic_terms, constant)
end

function MOI.get(
    model::Optimizer, ::MOI.ConstraintName,
    c::MOI.ConstraintIndex{MOI.ScalarQuadraticFunction{Float64}, S}
) where {S}
    return model[c].name
end

function MOI.set(
    model::Optimizer, ::MOI.ConstraintName,
    c::MOI.ConstraintIndex{MOI.ScalarQuadraticFunction{Float64}, S},
    name::String
) where {S}
    _update_if_necessary(model)
    set_strattrelement!(model.inner, "QCName", model[c].row, name)
    _require_update(model)
    model[c].name = name
    if model.name_to_constraint_index !== nothing && !haskey(model.name_to_constraint_index, name)
        model.name_to_constraint_index[c] = name
    else
        model.name_to_constraint_index = nothing
    end
    return
end

###
### VectorOfVariables-in-SOS{I|II}
###

const SOS = Union{MOI.SOS1{Float64}, MOI.SOS2{Float64}}

function Base.getindex(
    model::Optimizer,
    key::MOI.ConstraintIndex{MOI.VectorOfVariables, <:SOS}
)
    if haskey(model.sos_constraint_info, key.value)
        return model.sos_constraint_info[key.value]
    end
    throw(MOI.InvalidIndex(key))
end

sos_type(::MOI.SOS1) = :SOS1
sos_type(::MOI.SOS2) = :SOS2

function MOI.is_valid(
    model::Optimizer,
    c::MOI.ConstraintIndex{MOI.VectorOfVariables, S}
) where {S}
    info = get(model.sos_constraint_info, c.value, nothing)
    if info === nothing
        return false
    else
        return typeof(info.set) == S
    end
end

function MOI.add_constraint(
    model::Optimizer, f::MOI.VectorOfVariables, s::SOS
)
    columns = Int[model[v].column for v in f.variables]
    add_sos!(model.inner, sos_type(s), columns, s.weights)
    model.last_constraint_index += 1
    index = MOI.ConstraintIndex{MOI.VectorOfVariables, typeof(s)}(model.last_constraint_index)
    model.sos_constraint_info[index.value] = ConstraintInfo(
        length(model.sos_constraint_info) + 1, s
    )
    _require_update(model)
    return index
end

function MOI.delete(
    model::Optimizer, c::MOI.ConstraintIndex{MOI.VectorOfVariables, <:SOS}
)
    row = model[c].row
    _update_if_necessary(model)
    del_sos!(model.inner, [Cint(row)])
    _require_update(model)
    for (key, info) in model.sos_constraint_info
        if info.row > row
            info.row -= 1
        end
    end
    delete!(model.sos_constraint_info, c.value)
    return
end

function MOI.get(
    model::Optimizer, ::MOI.ConstraintName,
    c::MOI.ConstraintIndex{MOI.VectorOfVariables, <:Any}
)
    return model[c].name
end

function MOI.set(
    model::Optimizer, ::MOI.ConstraintName,
    c::MOI.ConstraintIndex{MOI.VectorOfVariables, <:Any}, name::String
)
    model[c].name = name
    if model.name_to_constraint_index !== nothing && !haskey(model.name_to_constraint_index, name)
        model.name_to_constraint_index[c] = name
    else
        model.name_to_constraint_index = nothing
    end
    return
end

function MOI.get(
    model::Optimizer, ::MOI.ConstraintSet,
    c::MOI.ConstraintIndex{MOI.VectorOfVariables, S}
) where {S <: SOS}
    _update_if_necessary(model)
    sparse_a, _ = get_sos(model.inner, model[c].row, 1)
    return S(sparse_a.nzval)
end

function MOI.get(
    model::Optimizer, ::MOI.ConstraintFunction,
    c::MOI.ConstraintIndex{MOI.VectorOfVariables, S}
) where {S <: SOS}
    _update_if_necessary(model)
    sparse_a, _ = get_sos(model.inner, model[c].row, 1)
    return MOI.VectorOfVariables(model.columns[sparse_a[1, :].nzind])
end

###
### Optimize methods.
###

function MOI.optimize!(model::Optimizer)
    # Note: Gurobi will call update regardless, so we don't have to.
    optimize(model.inner)
    model.has_infeasibility_cert =
    MOI.get(model, MOI.DualStatus()) == MOI.INFEASIBILITY_CERTIFICATE
    model.has_unbounded_ray =
        MOI.get(model, MOI.PrimalStatus()) == MOI.INFEASIBILITY_CERTIFICATE
    return
end

function MOI.get(model::Optimizer, ::MOI.TerminationStatus)
    stat = get_status(model.inner)
    if stat == :loaded
        return MOI.OPTIMIZE_NOT_CALLED
    elseif stat == :optimal
        return MOI.OPTIMAL
    elseif stat == :infeasible
        return MOI.INFEASIBLE
    elseif stat == :inf_or_unbd
        return MOI.INFEASIBLE_OR_UNBOUNDED
    elseif stat == :unbounded
        return MOI.DUAL_INFEASIBLE
    elseif stat == :cutoff
        return MOI.OBJECTIVE_LIMIT
    elseif stat == :iteration_limit
        return MOI.ITERATION_LIMIT
    elseif stat == :node_limit
        return MOI.NODE_LIMIT
    elseif stat == :time_limit
        return MOI.TIME_LIMIT
    elseif stat == :solution_limit
        return MOI.SOLUTION_LIMIT
    elseif stat == :interrupted
        return MOI.INTERRUPTED
    elseif stat == :numeric
        return MOI.NUMERICAL_ERROR
    elseif stat == :suboptimal
        return MOI.OTHER_LIMIT
    elseif stat == :inprogress
        return MOI.OTHER_ERROR
    elseif stat == :user_obj_limit
        return MOI.OBJECTIVE_LIMIT
    end
    return MOI.OTHER_ERROR
end

function MOI.get(model::Optimizer, ::MOI.PrimalStatus)
    stat = get_status(model.inner)
    if stat == :optimal
        return MOI.FEASIBLE_POINT
    elseif stat == :solution_limit
        return MOI.FEASIBLE_POINT
    elseif (stat == :inf_or_unbd || stat == :unbounded) && has_primal_ray(model)
        return MOI.INFEASIBILITY_CERTIFICATE
    elseif stat == :suboptimal
        return MOI.FEASIBLE_POINT
    elseif is_mip(model.inner) && get_sol_count(model.inner) > 0
        return MOI.FEASIBLE_POINT
    end
    return MOI.NO_SOLUTION
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

function MOI.get(model::Optimizer, ::MOI.DualStatus)
    stat = get_status(model.inner)
    if is_mip(model.inner)
        return MOI.NO_SOLUTION
    elseif stat == :optimal
        return MOI.FEASIBLE_POINT
    elseif stat == :solution_limit
        return MOI.FEASIBLE_POINT
    elseif (stat == :inf_or_unbd || stat == :infeasible) && has_dual_ray(model)
        return MOI.INFEASIBILITY_CERTIFICATE
    elseif stat == :suboptimal
        return MOI.FEASIBLE_POINT
    end
    return MOI.NO_SOLUTION
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

function MOI.get(model::Optimizer, ::MOI.VariablePrimal, x::MOI.VariableIndex)
    if model.has_unbounded_ray
        return get_dblattrelement(model.inner, "UnbdRay", model[x].column)
    else
        return get_dblattrelement(model.inner, "X", model[x].column)
    end
end

function MOI.get(
    model::Optimizer, ::MOI.ConstraintPrimal,
    c::MOI.ConstraintIndex{MOI.SingleVariable, <:Any}
)
    return MOI.get(model, MOI.VariablePrimal(), MOI.VariableIndex(c.value))
end

function MOI.get(
    model::Optimizer, ::MOI.ConstraintPrimal,
    c::MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64}, <:Any}
)
    row = model[c].row
    rhs = get_dblattrelement(model.inner, "RHS", row)
    slack = get_dblattrelement(model.inner, "Slack", row)
    return rhs - slack
end

function MOI.get(
    model::Optimizer, ::MOI.ConstraintPrimal,
    c::MOI.ConstraintIndex{MOI.ScalarQuadraticFunction{Float64}, <:Any}
)
    row = model[c].row
    rhs = get_dblattrelement(model.inner, "QCRHS", row)
    slack = get_dblattrelement(model.inner, "QCSlack", row)
    return rhs - slack
end

function dual_multiplier(model::Optimizer)
    return MOI.get(model, MOI.ObjectiveSense()) == MOI.MIN_SENSE ? 1.0 : -1.0
end

function MOI.get(
    model::Optimizer, ::MOI.ConstraintDual,
    c::MOI.ConstraintIndex{MOI.SingleVariable, MOI.LessThan{Float64}}
)
    x = get_dblattrelement(model.inner, "X", model[c].column)
    ub = get_dblattrelement(model.inner, "UB", model[c].column)
    if x ≈ ub
        return dual_multiplier(model) * get_dblattrelement(model.inner, "RC", model[c].column)
    else
        return 0.0
    end
end

function MOI.get(
    model::Optimizer, ::MOI.ConstraintDual,
    c::MOI.ConstraintIndex{MOI.SingleVariable, MOI.GreaterThan{Float64}}
)
    x = get_dblattrelement(model.inner, "X", model[c].column)
    lb = get_dblattrelement(model.inner, "LB", model[c].column)
    if x ≈ lb
        return dual_multiplier(model) * get_dblattrelement(model.inner, "RC", model[c].column)
    else
        return 0.0
    end
end

function MOI.get(
    model::Optimizer, ::MOI.ConstraintDual,
    c::MOI.ConstraintIndex{MOI.SingleVariable, MOI.EqualTo{Float64}}
)
    return dual_multiplier(model) * get_dblattrelement(model.inner, "RC", model[c].column)
end

function MOI.get(
    model::Optimizer, ::MOI.ConstraintDual,
    c::MOI.ConstraintIndex{MOI.SingleVariable, MOI.Interval{Float64}}
)
    return dual_multiplier(model) * get_dblattrelement(model.inner, "RC", model[c].column)
end

function MOI.get(
    model::Optimizer, ::MOI.ConstraintDual,
    c::MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64}, <:Any}
)
    # attr = model.has_infeasibility_cert ? "FarkasDual" : "Pi"
    # return dual_multiplier(model) * get_dblattrelement(model.inner, attr, model[c].row)
    if model.has_infeasibility_cert
        return dual_multiplier(model) * get_dblattrelement(model.inner, "FarkasDual", model[c].row)
    end
    return dual_multiplier(model) * get_dblattrelement(model.inner, "Pi", model[c].row)
end

function MOI.get(
    model::Optimizer, ::MOI.ConstraintDual,
    c::MOI.ConstraintIndex{MOI.ScalarQuadraticFunction{Float64}, <:Any}
)
    return dual_multiplier(model) * get_dblattrelement(model.inner, "QCPi", model[c].row)
end

MOI.get(model::Optimizer, ::MOI.ObjectiveValue) = get_dblattr(model.inner, "ObjVal")
MOI.get(model::Optimizer, ::MOI.ObjectiveBound) = get_dblattr(model.inner, "ObjBound")
MOI.get(model::Optimizer, ::MOI.SolveTime) = get_dblattr(model.inner, "RunTime")
MOI.get(model::Optimizer, ::MOI.SimplexIterations) = get_intattr(model.inner, "IterCount")
MOI.get(model::Optimizer, ::MOI.BarrierIterations) = get_intattr(model.inner, "BarIterCount")
MOI.get(model::Optimizer, ::MOI.NodeCount) = get_intattr(model.inner, "NodeCount")
MOI.get(model::Optimizer, ::MOI.RelativeGap) = get_dblattr(model.inner, "MIPGap")

MOI.supports(model::Optimizer, ::MOI.DualObjectiveValue) = true
MOI.get(model::Optimizer, ::MOI.DualObjectiveValue) = get_dblattr(model.inner, "ObjBound")

function MOI.get(model::Optimizer, ::MOI.ResultCount)
    if model.has_infeasibility_cert || model.has_unbounded_ray
        return 1
    end
    return get_intattr(model.inner, "SolCount")
end

function MOI.get(model::Optimizer, ::MOI.Silent)
    return model.silent
end

function MOI.set(model::Optimizer, ::MOI.Silent, flag::Bool)
    model.silent = flag
    output_flag = flag ? 0 : get(model.params, "OutputFlag", 1)
    setparam!(model.inner, "OutputFlag", output_flag)
    return
end

function MOI.get(model::Optimizer, ::MOI.Name)
    _update_if_necessary(model)
    return get_strattr(model.inner, "ModelName")
end

function MOI.set(model::Optimizer, ::MOI.Name, name::String)
    set_strattr!(model.inner, "ModelName", name)
    _require_update(model)
    return
end

MOI.get(model::Optimizer, ::MOI.NumberOfVariables) = length(model.variable_info)
function MOI.get(model::Optimizer, ::MOI.ListOfVariableIndices)
    return sort!(collect(keys(model.variable_info)), by = x -> x.value)
end

MOI.get(model::Optimizer, ::MOI.RawSolver) = model.inner

function MOI.set(
    model::Optimizer, ::MOI.VariablePrimalStart, x::MOI.VariableIndex, value::Float64
)
    Gurobi.set_dblattrelement!(model.inner, "Start", model[x].column, value)
    return
end

function MOI.get(model::Optimizer, ::MOI.VariablePrimalStart, x::MOI.VariableIndex)
    return get_dblattrelement(model.inner, "Start", model[x].column)
end

MOI.supports(::Optimizer, ::MOI.ConstraintPrimalStart) = false
MOI.supports(::Optimizer, ::MOI.ConstraintDualStart) = false

function MOI.get(model::Optimizer, ::MOI.NumberOfConstraints{F, S}) where {F, S}
    # TODO: this could be more efficient.
    return length(MOI.get(model, MOI.ListOfConstraintIndices{F, S}()))
end

bound_enums(::Type{<:MOI.LessThan}) = (LESS_THAN, LESS_AND_GREATER_THAN)
bound_enums(::Type{<:MOI.GreaterThan}) = (GREATER_THAN, LESS_AND_GREATER_THAN)
bound_enums(::Type{<:MOI.Interval}) = (INTERVAL,)
bound_enums(::Type{<:MOI.EqualTo}) = (EQUAL_TO,)
bound_enums(::Any) = (nothing,)

type_enums(::Type{MOI.ZeroOne}) = (BINARY,)
type_enums(::Type{MOI.Integer}) = (INTEGER,)
type_enums(::Type{<:MOI.Semicontinuous}) = (SEMICONTINUOUS,)
type_enums(::Type{<:MOI.Semiinteger}) = (SEMIINTEGER,)
type_enums(::Any) = (nothing,)

function MOI.get(
    model::Optimizer, ::MOI.ListOfConstraintIndices{MOI.SingleVariable, S}
) where {S}
    indices = MOI.ConstraintIndex{MOI.SingleVariable, S}[]
    for (key, info) in model.variable_info
        if info.bound in bound_enums(S) || info.type in type_enums(S)
            push!(indices, MOI.ConstraintIndex{MOI.SingleVariable, S}(key.value))
        end
    end
    return sort!(indices, by = x -> x.value)
end

function MOI.get(
    model::Optimizer,
    ::MOI.ListOfConstraintIndices{MOI.ScalarAffineFunction{Float64}, S}
) where {S}
    indices = MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64}, S}[]
    for (key, info) in model.affine_constraint_info
        if typeof(info.set) == S
            push!(indices, MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64}, S}(key))
        end
    end
    return sort!(indices, by = x -> x.value)
end

function MOI.get(
    model::Optimizer,
    ::MOI.ListOfConstraintIndices{MOI.ScalarQuadraticFunction{Float64}, S}
) where {S}
    indices = MOI.ConstraintIndex{MOI.ScalarQuadraticFunction{Float64}, S}[]
    for (key, info) in model.quadratic_constraint_info
        if typeof(info.set) == S
            push!(indices, MOI.ConstraintIndex{MOI.ScalarQuadraticFunction{Float64}, S}(key))
        end
    end
    return sort!(indices, by = x -> x.value)
end

function MOI.get(
    model::Optimizer, ::MOI.ListOfConstraintIndices{MOI.VectorOfVariables, S}
) where {S}
    indices = MOI.ConstraintIndex{MOI.VectorOfVariables, S}[]
    for (key, info) in model.sos_constraint_info
        if typeof(info.set) == S
            push!(indices, MOI.ConstraintIndex{MOI.VectorOfVariables, S}(key))
        end
    end
    return sort!(indices, by = x -> x.value)
end

function MOI.get(model::Optimizer, ::MOI.ListOfConstraints)
    constraints = Set{Any}()
    for info in values(model.variable_info)
        if info.bound == NONE
        elseif info.bound == LESS_THAN
            push!(constraints, (MOI.SingleVariable, MOI.LessThan{Float64}))
        elseif info.bound == GREATER_THAN
            push!(constraints, (MOI.SingleVariable, MOI.GreaterThan{Float64}))
        elseif info.bound == LESS_AND_GREATER_THAN
            push!(constraints, (MOI.SingleVariable, MOI.LessTHan{Float64}))
            push!(constraints, (MOI.SingleVariable, MOI.GreaterThan{Float64}))
        elseif info.bound == EQUAL_TO
            push!(constraints, (MOI.SingleVariable, MOI.EqualTo{Float64}))
        elseif info.bound == INTERVAL
            push!(constraints, (MOI.SingleVariable, MOI.Interval{Float64}))
        end
        if info.type == CONTINUOUS
        elseif info.type == BINARY
            push!(constraints, (MOI.SingleVariable, MOI.ZeroOne))
        elseif info.type == INTEGER
            push!(constraints, (MOI.SingleVariable, MOI.Integer))
        elseif info.type == SEMICONTINUOUS
            push!(constraints, (MOI.SingleVariable, MOI.Semicontinuous{Float64}))
        elseif info.type == SEMIINTEGER
            push!(constraints, (MOI.SingleVariable, MOI.Semiinteger{Float64}))
        end
    end
    for info in values(model.affine_constraint_info)
        push!(constraints, (MOI.ScalarAffineFunction{Float64}, typeof(info.set)))
    end
    for info in values(model.sos_constraint_info)
        push!(constraints, (MOI.VectorOfVariables, typeof(info.set)))
    end
    return collect(constraints)
end

function MOI.get(model::Optimizer, ::MOI.ObjectiveFunctionType)
    if model.objective_type == SCALAR_AFFINE
        return MOI.ScalarAffineFunction{Float64}
    elseif model.objective_type == SCALAR_QUADRATIC
        return MOI.ScalarQuadraticFunction{Float64}
    end
    error("Unknown objective type")
end

function MOI.modify(
    model::Optimizer,
    c::MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64}, <:Any},
    chg::MOI.ScalarCoefficientChange{Float64}
)
    chg_coeffs!(
        model.inner, model[c].row, model[chg.variable].column, chg.new_coefficient
    )
    _require_update(model)
end

function MOI.modify(
    model::Optimizer,
    c::MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}},
    chg::MOI.ScalarCoefficientChange{Float64}
)
    set_dblattrelement!(
        model.inner, "Obj", model[chg.variable].column, chg.new_coefficient
    )
    _require_update(model)
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
    replacement::MOI.ScalarAffineFunction, row::Int
)
    rows = fill(Cint(row), length(replacement.terms))
    cols = [Cint(model[t.variable_index].column) for t in replacement.terms]
    coefs = MOI.coefficient.(replacement.terms)
    chg_coeffs!(model.inner, rows, cols, coefs)
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
    replacement::MOI.ScalarAffineFunction, row::Int
)
    # First, zero out the old constraint function terms.
    rows = fill(Cint(row), length(previous.terms))
    cols = [Cint(model[t.variable_index].column) for t in previous.terms]
    coefs = fill(0.0, length(previous.terms))
    chg_coeffs!(model.inner, rows, cols, coefs)
    # Next, set the new constraint function terms.
    rows = fill(Cint(row), length(replacement.terms))
    cols = [Cint(model[t.variable_index].column) for t in replacement.terms]
    coefs = MOI.coefficient.(replacement.terms)
    chg_coeffs!(model.inner, rows, cols, coefs)
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
    f1::MOI.ScalarAffineFunction{Float64}, f2::MOI.ScalarAffineFunction{Float64}
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
    model::Optimizer, ::MOI.ConstraintFunction,
    c::MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64}, <:SCALAR_SETS},
    f::MOI.ScalarAffineFunction{Float64}
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
    row = model[c].row
    if _matching_sparsity_pattern(previous, replacement)
        _replace_with_matching_sparsity!(model, previous, replacement, row)
    else
        _replace_with_different_sparsity!(model, previous, replacement, row)
    end
    current_rhs = get_dblattrelement(model.inner, "RHS", row)
    new_rhs = current_rhs - (replacement.constant - previous.constant)
    set_dblattrelement!(model.inner, "RHS", row, new_rhs)
    _require_update(model)
    return
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

struct CallbackVariablePrimal <: MOI.AbstractVariableAttribute end

function load_callback_variable_primal(model, cb_data, cb_where)
    if cb_where != CB_MIPSOL
        error("`load_callback_variable_primal` must be called from `CB_MIPSOL`.")
    end
    resize!(model.callback_variable_primal, length(model.variable_info))
    cbget_mipsol_sol(cb_data, cb_where, model.callback_variable_primal)
    return
end

# Note: you must call load_callback_variable_primal first.
function MOI.get(
    model::Optimizer, ::CallbackVariablePrimal, x::MOI.VariableIndex
)
    return model.callback_variable_primal[model[x].column]
end

"""
    function cblazy!(
        cb_data::CallbackData, model::Optimizer,
        f::MOI.ScalarAffineFunction{Float64},
        s::Union{MOI.LessThan{Float64}, MOI.GreaterThan{Float64}, MOI.EqualTo{Float64}}
    )

Add a lazy cut to the model `m`.

You must have the option `LazyConstraints` set  via `Optimizer(LazyConstraint=1)`.
This can only be called in a callback from `CB_MIPSOL`.
"""
function cblazy!(
    cb_data::CallbackData, model::Optimizer,
    f::MOI.ScalarAffineFunction{Float64},
    s::Union{MOI.LessThan{Float64}, MOI.GreaterThan{Float64}, MOI.EqualTo{Float64}}
)
    indices, coefficients = indices_and_coefficients(model, f)
    sense, rhs = sense_and_rhs(s)
    return cblazy(cb_data, Cint.(indices), coefficients, Char(sense), rhs)
end

"""
    compute_conflict(model::Optimizer)

Compute a minimal subset of the constraints and variables that keep the model
infeasible.

See also `Gurobi.ConflictStatus` and `Gurobi.ConstraintConflictStatus`.

Note that if `model` is modified after a call to `compute_conflict`, the
conflict is not purged, and any calls to the above attributes will return values
for the original conflict without a warning.
"""
function compute_conflict(model::Optimizer)
    try
        computeIIS(model.inner)
    catch exc
        if isa(exc, GurobiError) && exc.code == 10015
            model.inner.conflict = Gurobi.GRB_INFEASIBLE
        else
            rethrow(exc)
        end
    end
    return
end

function _ensure_conflict_computed(model::Optimizer)
    if model.inner.conflict == -1
        error("Cannot access conflict status. Call `Gurobi.compute_conflict(model)` first. " *
              "In case the model is modified, the computed conflict will not be purged.")
    end
end

function _is_feasible(model::Optimizer)
    return model.inner.conflict == Gurobi.GRB_INFEASIBLE
end

"""
    ConflictStatus()

Return an `MOI.TerminationStatusCode` indicating the status of the last computed conflict.
If a minimal conflict is found, it will return `MOI.OPTIMAL`. If the problem is feasible, it will
return `MOI.INFEASIBLE`. If `compute_conflict` has not been called yet, it will return
`MOI.OPTIMIZE_NOT_CALLED`.
"""
struct ConflictStatus <: MOI.AbstractModelAttribute  end
MOI.is_set_by_optimize(::ConflictStatus) = true

function MOI.get(model::Optimizer, ::ConflictStatus)
    if model.inner.conflict == -1
        return MOI.OPTIMIZE_NOT_CALLED
    elseif model.inner.conflict == 0
        return MOI.OPTIMAL
    elseif model.inner.conflict == Gurobi.GRB_LOADED
        return MOI.OTHER_ERROR
    elseif model.inner.conflict == Gurobi.GRB_OPTIMAL
        return MOI.OPTIMAL
    elseif model.inner.conflict == Gurobi.GRB_INFEASIBLE
        return MOI.INFEASIBLE
    elseif model.inner.conflict == Gurobi.GRB_INF_OR_UNBD
        return MOI.INFEASIBLE_OR_UNBOUNDED
    elseif model.inner.conflict == Gurobi.GRB_USER_OBJ_LIMIT
        return MOI.OBJECTIVE_LIMIT
    elseif model.inner.conflict == Gurobi.GRB_ITERATION_LIMIT
        return MOI.ITERATION_LIMIT
    elseif model.inner.conflict == Gurobi.GRB_NODE_LIMIT
        return MOI.NODE_LIMIT
    elseif model.inner.conflict == Gurobi.GRB_TIME_LIMIT
        return MOI.TIME_LIMIT
    elseif model.inner.conflict == Gurobi.GRB_SOLUTION_LIMIT
        return MOI.SOLUTION_LIMIT
    elseif model.inner.conflict == Gurobi.GRB_INTERRUPTED
        return MOI.INTERRUPTED
    elseif model.inner.conflict == Gurobi.GRB_NUMERIC
        return MOI.NUMERICAL_ERROR
    elseif model.inner.conflict == Gurobi.GRB_SUBOPTIMAL
        return MOI.OTHER_LIMIT
    elseif model.inner.conflict == Gurobi.GRB_INPROGRESS
        return MOI.OTHER_ERROR
    else
        return MOI.OTHER_ERROR
    end
end

function MOI.supports(::Optimizer, ::ConflictStatus)
    return true
end

"""
    ConstraintConflictStatus()
A Boolean constraint attribute indicating whether the constraint participates in the last computed conflict.
"""
struct ConstraintConflictStatus <: MOI.AbstractConstraintAttribute end
MOI.is_set_by_optimize(::ConstraintConflictStatus) = true

function MOI.get(model::Optimizer, ::ConstraintConflictStatus, index::MOI.ConstraintIndex{<:MOI.SingleVariable, <:LQOI.LE})
    _ensure_conflict_computed(model)
    return !_is_feasible(model) && Bool(get_intattrelement(model.inner, "IISUB", LQOI.get_column(model, model[index])))
end

function MOI.get(model::Optimizer, ::ConstraintConflictStatus, index::MOI.ConstraintIndex{<:MOI.SingleVariable, <:LQOI.GE})
    _ensure_conflict_computed(model)
    return !_is_feasible(model) && Bool(get_intattrelement(model.inner, "IISLB", LQOI.get_column(model, model[index])))
end

function MOI.get(model::Optimizer, ::ConstraintConflictStatus, index::MOI.ConstraintIndex{<:MOI.SingleVariable, <:Union{LQOI.EQ, LQOI.IV}})
    _ensure_conflict_computed(model)
    return !_is_feasible(model) && (
        Bool(get_intattrelement(model.inner, "IISUB", LQOI.get_column(model, model[index]))) || Bool(get_intattrelement(model.inner, "IISLB", model[index])))
end

function MOI.get(model::Optimizer, ::ConstraintConflictStatus, index::MOI.ConstraintIndex{<:MOI.ScalarAffineFunction, <:Union{LQOI.LE, LQOI.GE, LQOI.EQ}})
    _ensure_conflict_computed(model)
    return !_is_feasible(model) && Bool(get_intattrelement(model.inner, "IISConstr", model[index]))
end

function MOI.get(model::Optimizer, ::ConstraintConflictStatus, index::MOI.ConstraintIndex{<:MOI.ScalarQuadraticFunction, <:Union{LQOI.LE, LQOI.GE}})
    _ensure_conflict_computed(model)
    return !_is_feasible(model) && Bool(get_intattrelement(model.inner, "IISQConstr", model[index]))
end

function MOI.supports(::Optimizer, ::ConstraintConflictStatus, ::Type{MOI.ConstraintIndex{<:MOI.SingleVariable, <:LQOI.LinSets}})
    return true
end

function MOI.supports(::Optimizer, ::ConstraintConflictStatus, ::Type{MOI.ConstraintIndex{<:MOI.ScalarAffineFunction, <:Union{LQOI.LE, LQOI.GE, LQOI.EQ}}})
    return true
end

function MOI.supports(::Optimizer, ::ConstraintConflictStatus, ::Type{MOI.ConstraintIndex{<:MOI.ScalarQuadraticFunction, <:Union{LQOI.LE, LQOI.GE}}})
    return true
end
