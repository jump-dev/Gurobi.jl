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

function _add_expression_tree_node(
    model::Optimizer,
    stack::Vector{Tuple{Any,Cint}},
    opcode::Vector{Cint},
    data::Vector{Cdouble},
    parent::Vector{Cint},
    s::MOI.ScalarNonlinearFunction,
    parent_index::Cint,
)
    # Append an operator to opcode/data/parent arrays and push new children onto
    # the stack
    current_opcode = get(_OPCODE_MAP, s.head, nothing)
    if current_opcode === nothing
        throw(MOI.UnsupportedNonlinearOperator(s.head))
    elseif current_opcode == GRB_OPCODE_MINUS
        # MINUS opcode is binary only
        if length(s.args) == 1
            current_opcode = GRB_OPCODE_UMINUS
        else
            @assert length(s.args) == 2
        end
    end
    append!(opcode, current_opcode)
    append!(data, -1.0)
    append!(parent, parent_index)
    # Children all have the same parent (the entry just added)
    for expr in reverse(s.args)
        push!(stack, (expr, length(opcode) - 1))
    end
    return
end

function _add_expression_tree_node(
    model::Optimizer,
    stack::Vector{Tuple{Any,Cint}},
    opcode::Vector{Cint},
    data::Vector{Cdouble},
    parent::Vector{Cint},
    s::MOI.VariableIndex,
    parent_index::Cint,
)
    # Variable leaf node: just append to the arrays, nothing new goes on the stack.
    append!(opcode, GRB_OPCODE_VARIABLE)
    append!(data, c_column(model, s))
    append!(parent, parent_index)
    return
end

function _add_expression_tree_node(
    model::Optimizer,
    stack::Vector{Tuple{Any,Cint}},
    opcode::Vector{Cint},
    data::Vector{Cdouble},
    parent::Vector{Cint},
    s::Union{Float64,Int},
    parent_index::Cint,
)
    # Constant leaf node: just append to the arrays, nothing new goes on the stack.
    append!(opcode, GRB_OPCODE_CONSTANT)
    append!(data, s)
    append!(parent, parent_index)
    return
end

function _add_expression_tree_node(
    model::Optimizer,
    stack::Vector{Tuple{Any,Cint}},
    opcode::Vector{Cint},
    data::Vector{Cdouble},
    parent::Vector{Cint},
    term::MOI.ScalarAffineTerm{Float64},
    parent_index::Cint,
)
    # Add a binary MULTIPLY node with a constant and a variable as children
    append!(opcode, GRB_OPCODE_MULTIPLY)
    append!(data, -1.0)
    append!(parent, parent_index)
    multiply_parent_index = Cint(length(opcode) - 1)
    _add_expression_tree_node(
        model,
        stack,
        opcode,
        data,
        parent,
        term.coefficient,
        multiply_parent_index,
    )
    _add_expression_tree_node(
        model,
        stack,
        opcode,
        data,
        parent,
        term.variable,
        multiply_parent_index,
    )
    return
end

function _add_expression_tree_node(
    model::Optimizer,
    stack::Vector{Tuple{Any,Cint}},
    opcode::Vector{Cint},
    data::Vector{Cdouble},
    parent::Vector{Cint},
    s::MOI.ScalarAffineFunction{Float64},
    parent_index::Cint,
)
    # Add an N-ary PLUS node with all terms as children. Nothing new goes on the
    # stack, the expression is expanded in-place.
    append!(opcode, GRB_OPCODE_PLUS)
    append!(data, -1.0)
    append!(parent, parent_index)
    plus_parent_index = Cint(length(opcode) - 1)
    if !iszero(s.constant)
        _add_expression_tree_node(
            model,
            stack,
            opcode,
            data,
            parent,
            s.constant,
            plus_parent_index,
        )
    end
    for term in s.terms
        _add_expression_tree_node(
            model,
            stack,
            opcode,
            data,
            parent,
            term,
            plus_parent_index,
        )
    end
    return
end

function _add_expression_tree_node(
    model::Optimizer,
    stack::Vector{Tuple{Any,Cint}},
    opcode::Vector{Cint},
    data::Vector{Cdouble},
    parent::Vector{Cint},
    term::MOI.ScalarQuadraticTerm{Float64},
    parent_index::Cint,
)
    # Add a binary MULTIPLY node with a constant and two variables as children
    append!(opcode, GRB_OPCODE_MULTIPLY)
    append!(data, -1.0)
    append!(parent, parent_index)
    multiply_parent_index = Cint(length(opcode) - 1)
    # https://jump.dev/MathOptInterface.jl/stable/reference/standard_form
    # ScalarQuadraticFunction stores diagonal ScalarQuadraticTerms multiplied by
    # 2
    coeff = term.coefficient
    if term.variable_1 == term.variable_2
        coeff /= 2
    end
    if !isone(coeff)
        _add_expression_tree_node(
            model,
            stack,
            opcode,
            data,
            parent,
            coeff,
            multiply_parent_index,
        )
    end
    _add_expression_tree_node(
        model,
        stack,
        opcode,
        data,
        parent,
        term.variable_1,
        multiply_parent_index,
    )
    _add_expression_tree_node(
        model,
        stack,
        opcode,
        data,
        parent,
        term.variable_2,
        multiply_parent_index,
    )
    return
end

function _add_expression_tree_node(
    model::Optimizer,
    stack::Vector{Tuple{Any,Cint}},
    opcode::Vector{Cint},
    data::Vector{Cdouble},
    parent::Vector{Cint},
    s::MOI.ScalarQuadraticFunction{Float64},
    parent_index::Cint,
)
    # Add an N-ary PLUS node with all terms as children. Nothing new goes on the
    # stack, the expression is expanded in-place.
    append!(opcode, GRB_OPCODE_PLUS)
    append!(data, -1.0)
    append!(parent, parent_index)
    plus_parent_index = Cint(length(opcode) - 1)
    if !iszero(s.constant)
        _add_expression_tree_node(
            model,
            stack,
            opcode,
            data,
            parent,
            s.constant,
            plus_parent_index,
        )
    end
    for term in s.affine_terms
        _add_expression_tree_node(
            model,
            stack,
            opcode,
            data,
            parent,
            term,
            plus_parent_index,
        )
    end
    for term in s.quadratic_terms
        _add_expression_tree_node(
            model,
            stack,
            opcode,
            data,
            parent,
            term,
            plus_parent_index,
        )
    end
    return
end

function _process_nonlinear(
    model::Optimizer,
    f,
    opcode::Vector{Cint},
    data::Vector{Cdouble},
    parent::Vector{Cint},
)
    stack = Tuple{Any,Cint}[(f, Cint(-1))]
    while !isempty(stack)
        s, parent_index = pop!(stack)
        parent_index = _add_expression_tree_node(
            model,
            stack,
            opcode,
            data,
            parent,
            s,
            parent_index,
        )
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
    _process_nonlinear(model, f, opcode, data, parent)
    # Add resultant variable. We don't use MOI.add_constrained_variable because
    # we don't want it to show up in the bound constraints, etc.
    column = _get_next_column(model)
    lb, ub = _bounds(s)
    lb = something(lb, -Inf)
    ub = something(ub, Inf)
    ret = GRBaddvar(model, 0, C_NULL, C_NULL, 0.0, lb, ub, GRB_CONTINUOUS, "")
    _check_ret(model, ret)
    ret = GRBaddgenconstrNL(
        model,
        C_NULL,
        column - 1,
        length(opcode),
        opcode,
        data,
        parent,
    )
    _check_ret(model, ret)
    _require_update(model, model_change = true)
    model.last_constraint_index += 1
    model.nl_constraint_info[model.last_constraint_index] =
        _NLConstraintInfo(length(model.nl_constraint_info) + 1, s, column)
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
    c::MOI.ConstraintIndex{MOI.ScalarNonlinearFunction},
)
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
    # Delete resultant variable from the Gurobi model. These are not tracked in
    # model.variable_info but they do need to be accounted for in index
    # adjustment.
    del_cols = [Cint(info.resvar_index - 1)]
    ret = GRBdelvars(model, length(del_cols), del_cols)
    _check_ret(model, ret)
    append!(model.columns_deleted_since_last_update, del_cols .+ 1)
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
    # Delete resultant variables from the Gurobi model for all removed
    # constraints. These are not tracked in model.variable_info but they do
    # need to be accounted for in index adjustment.
    del_cols = [Cint(_info(model, c).resvar_index - 1) for c in cs]
    ret = GRBdelvars(model, length(del_cols), del_cols)
    _check_ret(model, ret)
    append!(model.columns_deleted_since_last_update, del_cols .+ 1)
    _require_update(model, model_change = true)
    # Remove entries from nl constraint tracking
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

function MOI.get(
    model::Optimizer,
    attr::MOI.ConstraintPrimal,
    c::MOI.ConstraintIndex{MOI.ScalarNonlinearFunction},
)
    _throw_if_optimize_in_progress(model, attr)
    MOI.check_result_index_bounds(model, attr)
    info = _info(model, c)
    key = "X"
    if attr.result_index > 1
        MOI.set(
            model,
            MOI.RawOptimizerAttribute("SolutionNumber"),
            attr.result_index - 1,
        )
        key = "Xn"
    end
    valueP = Ref{Cdouble}()
    ret = GRBgetdblattrelement(model, key, Cint(info.resvar_index - 1), valueP)
    _check_ret(model, ret)
    return valueP[]
    return
end
