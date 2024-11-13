# Copyright (c) 2015 Dahua Lin, Miles Lubin, Joey Huchette, Iain Dunning, and
#   contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

_supports_nonlinear() = _GUROBI_VERSION >= v"12.0.0"

if _supports_nonlinear()
    const _OPCODE_MAP = Dict(
        :+ => GRB_OPCODE_PLUS,
        :- => GRB_OPCODE_MINUS,
        :* => GRB_OPCODE_MULTIPLY,
        :/ => GRB_OPCODE_DIVIDE,
        # GRB_OPCODE_UMINUS        6
        # GRB_OPCODE_SQUARE        7
        :sqrt => GRB_OPCODE_SQRT,
        :sin => GRB_OPCODE_SIN,
        :cos => GRB_OPCODE_COS,
        :tan => GRB_OPCODE_TAN,
        :^ => GRB_OPCODE_POW,
        :exp => GRB_OPCODE_EXP,
        :log => GRB_OPCODE_LOG,
        :log2 => GRB_OPCODE_LOG2,
        :log10 => GRB_OPCODE_LOG10,
        :logistic => GRB_OPCODE_LOGISTIC,
    )
else
    const _OPCODE_MAP = Dict()
end

function MOI.supports_constraint(
    ::Optimizer,
    ::Type{MOI.ScalarNonlinearFunction},
    ::Type{
        <:Union{
            MOI.LessThan{Float64},
            MOI.GreaterThan{Float64},
            MOI.EqualTo{Float64},
        },
    },
)
    return _supports_nonlinear()
end

function _info(
    model::Optimizer,
    c::MOI.ConstraintIndex{MOI.ScalarNonlinearFunction},
)
    if haskey(model.nl_constraint_info, c.value)
        return model.nl_constraint_info[c.value]
    end
    return throw(MOI.InvalidIndex(c))
end

function _add_expression_tree_constant(
    opcode::Vector{Cint},
    data::Vector{Cdouble},
    parent::Vector{Cint},
    coeff::Float64,
    index::Cint,
)
    append!(opcode, GRB_OPCODE_CONSTANT)
    append!(data, coeff)
    append!(parent, index)
    return
end

function _add_expression_tree_variable(
    opcode::Vector{Cint},
    data::Vector{Cdouble},
    parent::Vector{Cint},
    coeff::Float64,
    var_index::Cint,
    current_index::Cint,
    parent_index::Cint,
)
    if !isone(coeff)
        append!(opcode, GRB_OPCODE_MULTIPLY)
        append!(data, -1.0)
        append!(parent, parent_index)
        _add_expression_tree_constant(
            opcode,
            data,
            parent,
            coeff,
            current_index,
        )
        append!(opcode, GRB_OPCODE_VARIABLE)
        append!(data, var_index)
        append!(parent, current_index)
        current_index += Cint(2)
        return current_index
    else
        append!(data, var_index)
        append!(opcode, GRB_OPCODE_VARIABLE)
        append!(parent, parent_index)
        current_index += Cint(1)
        return current_index
    end
end

function _process_nonlinear(
    model::Optimizer,
    f,
    opcode::Vector{Cint},
    data::Vector{Cdouble},
    parent::Vector{Cint},
)
    stack = Vector{Tuple{Any,Cint}}([(f, Cint(-1))])
    current_index = Cint(-1)
    while !isempty(stack)
        current_index += Cint(1)
        s, parent_index = pop!(stack)
        if s isa MOI.ScalarNonlinearFunction
            ret = get(_OPCODE_MAP, s.head, nothing)
            if ret === nothing
                throw(MOI.UnsupportedNonlinearOperator(s.head))
            elseif s.head == :- && length(s.args) == 1  # Special handling for unary -
                append!(opcode, GRB_OPCODE_UMINUS)
                append!(data, -1.0)
                append!(parent, parent_index)
            else
                append!(opcode, ret)
                append!(data, -1.0)
                append!(parent, parent_index)
            end
            for expr in reverse(s.args)
                push!(stack, (expr, current_index))
            end
        elseif s isa MOI.VariableIndex
            _add_expression_tree_variable(
                opcode,
                data,
                parent,
                1.0,
                c_column(model, s),
                current_index,
                parent_index,
            )
        elseif s isa MOI.ScalarAffineFunction{Float64}
            if length(s.terms) > 1
                append!(opcode, GRB_OPCODE_PLUS)
                append!(data, -1.0)
                append!(parent, parent_index)
                parent_index += Cint(1)
            end
            if !iszero(s.constant)
                append!(opcode, GRB_OPCODE_CONSTANT)
                append!(data, s.constant)
                append!(parent, parent_index)
                current_index += Cint(1)
            end
            for term in s.terms
                var_index = c_column(model, term.variable)
                coeff = term.coefficient
                current_index = _add_expression_tree_variable(
                    opcode,
                    data,
                    parent,
                    coeff,
                    var_index,
                    current_index,
                    parent_index,
                )
                current_index += Cint(1)
            end
        elseif s isa Float64
            _add_expression_tree_constant(opcode, data, parent, s, parent_index)
        elseif s isa Int
            _add_expression_tree_constant(
                opcode,
                data,
                parent,
                convert(Float64, s),
                parent_index,
            )
        else
            throw(MOI.UnsupportedAttribute(s))
        end
    end
    return
end

function MOI.add_constraint(
    model::Optimizer,
    f::MOI.ScalarNonlinearFunction,
    s::_SCALAR_SETS,
)
    opcode = Cint[]
    data = Cdouble[]
    parent = Cint[]
    sense, rhs = _sense_and_rhs(s)
    _process_nonlinear(model, f, opcode, data, parent)
    # Add resultant variable
    vi, ci = MOI.add_constrained_variable(model, s)
    resvar_index = c_column(model, vi)
    ret = GRBaddgenconstrNL(
        model,
        C_NULL,
        resvar_index,
        length(opcode),
        opcode,
        data,
        parent,
    )
    _check_ret(model, ret)
    _require_update(model, model_change = true)
    model.last_constraint_index += 1
    model.nl_constraint_info[model.last_constraint_index] =
        _NLConstraintInfo(length(model.nl_constraint_info) + 1, s, vi)
    return MOI.ConstraintIndex{MOI.ScalarNonlinearFunction,typeof(s)}(
        model.last_constraint_index,
    )
end

function MOI.is_valid(
    model::Optimizer,
    c::MOI.ConstraintIndex{MOI.ScalarNonlinearFunction,S},
) where {S}
    info = get(model.nl_constraint_info, c.value, nothing)
    return info !== nothing && info.set isa S
end

function MOI.delete(
    model::Optimizer,
    c::MOI.ConstraintIndex{MOI.ScalarNonlinearFunction,S},
) where {S}
    info = _info(model, c)
    _update_if_necessary(model)
    ret = GRBdelgenconstrs(model, 1, Ref(Cint(info.row - 1)))
    _check_ret(model, ret)
    for (_, info_2) in model.nl_constraint_info
        if info_2.row > info.row
            info_2.row -= 1
        end
    end
    delete!(model.nl_constraint_info, c.value)
    model.name_to_constraint_index = nothing
    # Remove resultant variable
    MOI.delete(model, info.resvar)
    _require_update(model, model_change = true)
    return
end

function MOI.delete(
    model::Optimizer,
    cs::Vector{<:MOI.ConstraintIndex{MOI.ScalarNonlinearFunction}},
)
    rows_to_delete = sort!([Cint(_info(model, c).row - 1) for c in cs])
    _update_if_necessary(model)
    ret = GRBdelgenconstrs(model, length(rows_to_delete), rows_to_delete)
    _check_ret(model, ret)
    for (_, info) in model.nl_constraint_info
        info.row -= searchsortedlast(rows_to_delete, info.row - 1)
    end
    model.name_to_constraint_index = nothing
    # Delete resultant variables
    resvars = [_info(model, c).resvar for c in cs]
    MOI.delete(model, resvars)
    cs_values = sort!(getfield.(cs, :value))
    filter!(model.nl_constraint_info) do pair
        return isempty(searchsorted(cs_values, pair.first))
    end
    _require_update(model, model_change = true)
    return
end

function MOI.get(
    model::Optimizer,
    ::MOI.ConstraintName,
    c::MOI.ConstraintIndex{MOI.ScalarNonlinearFunction},
)
    return _info(model, c).name
end

function MOI.set(
    model::Optimizer,
    ::MOI.ConstraintName,
    c::MOI.ConstraintIndex{MOI.ScalarNonlinearFunction},
    name::String,
)
    _update_if_necessary(model, force = true)
    info = _info(model, c)
    info.name = name
    if length(name) <= GRB_MAX_NAMELEN
        ret = GRBsetstrattrelement(
            model,
            "GenConstrName",
            Cint(info.row - 1),
            name,
        )
        _check_ret(model, ret)
        _require_update(model, attribute_change = true)
    end
    model.name_to_constraint_index = nothing
    return
end
