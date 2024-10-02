# Copyright (c) 2015 Dahua Lin, Miles Lubin, Joey Huchette, Iain Dunning, and
#   contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

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
    return _GUROBI_VERSION >= v"12.0.0"
end

function _info(
    model::Optimizer,
    c::MOI.ConstraintIndex{MOI.ScalarNonlinearFunction,S},
) where {S}
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
    return append!(parent, index)
end

function _add_expression_tree_variable(
    opcode::Vector{Cint},
    data::Vector{Cdouble},
    parent::Vector{Cint},
    coeff::Float64,
    var::MOI.VariableIndex,
    var_index::Cint,
    current_index::Cint,
    parent_index::Cint,
)
    if coeff > 1.0
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
    else
        append!(data, var_index)
        append!(opcode, GRB_OPCODE_VARIABLE)
        append!(parent, parent_index)
    end
end

# Check if a nonlinear is actually just a constant
function _check_nonlinear_constant(expr)
    return (
        typeof(expr) == MOI.ScalarNonlinearFunction &&
        expr.head in [:+, :-] &&
        length(expr.args) == 1 &&
        typeof(expr.args[1]) == Float64
    )
end

# Check if a nonlinear is actually just an affine expression
# Cases:
#   1. constant * var
#   2. +/- var
#   3. +/- constant * var
function _check_nonlinear_affine(expr)
    return typeof(expr) == MOI.ScalarNonlinearFunction && (
        ( # Case 1.
            expr.head == :* &&
            length(expr.args) == 2 &&
            typeof(expr.args[2]) == MOI.VariableIndex
        ) || ( # Cases 2. and 3.
            expr.head in [:+, :-] &&
            length(expr.args) == 1 &&
            (
                typeof(expr.args[1]) == MOI.VariableIndex ||
                _check_nonlinear_affine(expr.args[1])
            )
        )
    )
end

function _process_nonlinear(
    model::Optimizer,
    f,
    opcode::Vector{Cint},
    data::Vector{Cdouble},
    parent::Vector{Cint},
)
    # TODO: use type hints here instead of Any
    stack = Vector{Tuple{Any,Cint}}([(f, Cint(-1))])
    current_index = Cint(-1)

    while length(stack) != 0
        current_index += Cint(1)
        s, parent_index = pop!(stack)

        if typeof(s) == MOI.ScalarNonlinearFunction
            ret = get(_OPCODE_MAP, s.head, nothing)
            if ret === nothing
                throw(MOI.UnsupportedNonlinearOperator(s.head))
            elseif s.head == :- && length(s.args) == 1  # Special handling for unary -
                append!(opcode, GRB_OPCODE_UMINUS)
            elseif !_check_nonlinear_affine(s) && !_check_nonlinear_constant(s)
                append!(opcode, ret)
                append!(data, -1.0)
                append!(parent, parent_index)
            end

            for expr in reverse(s.args)
                push!(stack, (expr, current_index))
            end
        elseif typeof(s) == MOI.VariableIndex
            _add_expression_tree_variable(
                opcode,
                data,
                parent,
                1.0,
                s,
                c_column(model, s),
                current_index,
                parent_index,
            )
        elseif typeof(s) == MOI.ScalarAffineFunction{Float64}
            for term in s.terms
                var_index = c_column(model, term.variable)
                coeff = term.coefficient
                _add_expression_tree_variable(
                    opcode,
                    data,
                    parent,
                    coeff,
                    term.variable,
                    var_index,
                    current_index,
                    parent_index,
                )
            end
        elseif typeof(s) == Float64
            _add_expression_tree_constant(opcode, data, parent, s, parent_index)
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
    opcode = Vector{Cint}()
    data = Vector{Cdouble}()
    parent = Vector{Cint}()

    f_with_rhs = f
    sense, rhs = _sense_and_rhs(s)

    if rhs > 0
        f_with_rhs = MOI.ScalarNonlinearFunction(:-, Any[f, rhs])
    elseif rhs < 0
        f_with_rhs = MOI.ScalarNonlinearFunction(:+, Any[f, rhs])
    end

    _process_nonlinear(model, f_with_rhs, opcode, data, parent)

    # Add resultant variable
    vi, ci = MOI.add_constrained_variable(model, s)

    ret = GRBaddgenconstrNL(
        model,
        C_NULL,
        c_column(model, vi),
        length(opcode),
        opcode,
        data,
        parent,
    )
    _check_ret(model, ret)
    _require_update(model, model_change = true)
    model.last_constraint_index += 1
    model.nl_constraint_info[model.last_constraint_index] =
        _ConstraintInfo(length(model.nl_constraint_info) + 1, s)
    return MOI.ConstraintIndex{MOI.ScalarNonlinearFunction,typeof(s)}(
        model.last_constraint_index,
    )
end

function MOI.is_valid(
    model::Optimizer,
    c::MOI.ConstraintIndex{MOI.ScalarNonlinearFunction,S},
) where {S}
    info = get(model.nl_constraint_info, c.value, nothing)
    return info !== nothing && typeof(info.set) == S
end

function MOI.delete(
    model::Optimizer,
    c::MOI.ConstraintIndex{MOI.ScalarNonlinearFunction,S},
) where {S}
    info = _info(model, c)
    ret = GRBdelgenconstrs(model, 1, Ref{Cint}(info.row - 1))
    _check_ret(model, ret)
    for (key, info_2) in model.nl_constraint_info
        if info_2.row > info.row
            info_2.row -= 1
        end
    end
    delete!(model.nl_constraint_info, c.value)
    model.name_to_constraint_index = nothing
    _update_if_necessary(model, force = true)
    return
end



# TODO: all the functions below

function MOI.get(
    model::Optimizer,
    ::MOI.ConstraintSet,
    c::MOI.ConstraintIndex{MOI.ScalarNonlinearFunction,S},
) where {S}
    throw("not implemented")
end

function MOI.get(
    model::Optimizer,
    ::MOI.ConstraintFunction,
    c::MOI.ConstraintIndex{MOI.ScalarNonlinearFunction,S},
) where {S}
    throw("not implemented")
end

function MOI.get(
    model::Optimizer,
    ::MOI.ConstraintName,
    c::MOI.ConstraintIndex{MOI.ScalarNonlinearFunction,S},
) where {S}
    return _info(model, c).name
end

function MOI.set(
    model::Optimizer,
    ::MOI.ConstraintName,
    c::MOI.ConstraintIndex{MOI.ScalarNonlinearFunction,S},
    name::String,
) where {S}
    throw("not implemented")
    # info = _info(model, c)
    # info.name = name
    # model.name_to_constraint_index = nothing
    return
end
