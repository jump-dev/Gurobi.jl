function _info(
    model::Optimizer,
    c::MOI.ConstraintIndex{MOI.VectorAffineFunction{Float64},<:MOI.Indicator},
)
    if haskey(model.indicator_constraint_info, c.value)
        return model.indicator_constraint_info[c.value]
    end
    return throw(MOI.InvalidIndex(c))
end

function MOI.supports_constraint(
    ::Optimizer,
    ::Type{MOI.VectorAffineFunction{Float64}},
    ::Type{<:MOI.Indicator{A,S}},
) where {A,S<:_SUPPORTED_SCALAR_SETS}
    return true
end

function MOI.is_valid(
    model::Optimizer,
    c::MOI.ConstraintIndex{MOI.VectorAffineFunction{Float64},S},
) where {S<:MOI.Indicator}
    info = get(model.indicator_constraint_info, c.value, nothing)
    if info === nothing
        return false
    end
    return isa(info.set, S)
end

function MOI.get(
    model::Optimizer,
    ::MOI.ConstraintSet,
    c::MOI.ConstraintIndex{<:MOI.VectorAffineFunction,<:MOI.Indicator},
)
    MOI.throw_if_not_valid(model, c)
    return _info(model, c).set
end

function MOI.get(
    model::Optimizer,
    ::MOI.ConstraintFunction,
    c::MOI.ConstraintIndex{<:MOI.VectorAffineFunction,<:MOI.Indicator},
)
    MOI.throw_if_not_valid(model, c)
    _update_if_necessary(model)
    info = _info(model, c)
    binvarP, nvarsP = Ref{Cint}(), Ref{Cint}()
    ret = GRBgetgenconstrIndicator(
        model,
        Cint(info.row - 1),
        binvarP,
        C_NULL,
        nvarsP,
        C_NULL,
        C_NULL,
        C_NULL,
        C_NULL,
    )
    _check_ret(model, ret)
    vars = Vector{Cint}(undef, nvarsP[])
    vals = Vector{Cdouble}(undef, nvarsP[])
    rhsP = Ref{Cdouble}()
    ret = GRBgetgenconstrIndicator(
        model,
        Cint(info.row - 1),
        C_NULL,
        C_NULL,
        C_NULL,
        vars,
        vals,
        C_NULL,
        rhsP,
    )
    _check_ret(model, ret)
    terms = Vector{MOI.VectorAffineTerm{Float64}}(undef, nvarsP[] + 1)
    x = model.variable_info[CleverDicts.LinearIndex(binvarP[] + 1)].index
    terms[1] = MOI.VectorAffineTerm(1, MOI.ScalarAffineTerm(1.0, x))
    for i in 1:nvarsP[]
        x = model.variable_info[CleverDicts.LinearIndex(vars[i] + 1)].index
        terms[i+1] = MOI.VectorAffineTerm(2, MOI.ScalarAffineTerm(vals[i], x))
    end
    _, rhs = _sense_and_rhs(info.set.set)
    return MOI.VectorAffineFunction(terms, [0.0, rhs - rhsP[]])
end

function MOI.add_constraint(
    model::Optimizer,
    func::MOI.VectorAffineFunction{Float64},
    s::MOI.Indicator{A,S},
) where {A,S<:_SUPPORTED_SCALAR_SETS}
    if !iszero(func.constants[1])
        error(
            "Constant in output_index 1 should be 0. Got $(func.constants[1])",
        )
    end
    vars = Vector{Cint}(undef, length(func.terms) - 1)
    vals = Vector{Cdouble}(undef, length(func.terms) - 1)
    binvar = Cint(-1)
    i = 1
    for vector_term in func.terms
        row = vector_term.output_index
        term = vector_term.scalar_term
        if row == 1
            if binvar !== Cint(-1)
                error("There should be exactly one term in output_index 1")
            elseif !isapprox(term.coefficient, 1.0)
                error(
                    "Expected coefficient in front of indicator variable to " *
                    "be 1.0. Got $(term.coefficient).",
                )
            end
            binvar = Cint(column(model, term.variable) - 1)
        else
            @assert row == 2
            vars[i] = Cint(column(model, term.variable) - 1)
            vals[i] = term.coefficient
            i += 1
        end
    end
    sense, rhs = _sense_and_rhs(s.set)
    rhs -= func.constants[2]
    ret = GRBaddgenconstrIndicator(
        model,
        "",
        binvar,
        A == MOI.ACTIVATE_ON_ONE ? Cint(1) : Cint(0),
        length(vars),
        vars,
        vals,
        sense,
        rhs,
    )
    _check_ret(model, ret)
    model.last_constraint_index += 1
    info = _ConstraintInfo(length(model.indicator_constraint_info) + 1, s)
    model.indicator_constraint_info[model.last_constraint_index] = info
    _require_update(model)
    return MOI.ConstraintIndex{typeof(func),typeof(s)}(
        model.last_constraint_index,
    )
end

function MOI.get(
    model::Optimizer,
    ::MOI.ListOfConstraintIndices{<:MOI.VectorAffineFunction,S},
) where {S<:MOI.Indicator}
    indices = MOI.ConstraintIndex{MOI.VectorAffineFunction{Float64},S}[
        MOI.ConstraintIndex{MOI.VectorAffineFunction{Float64},S}(key) for
        (key, info) in model.indicator_constraint_info if isa(info.set, S)
    ]
    return sort!(indices, by = x -> x.value)
end

function MOI.get(
    model::Optimizer,
    ::MOI.NumberOfConstraints{MOI.VectorAffineFunction{Float64},S},
) where {S<:MOI.Indicator}
    return count(x -> isa(x.set, S), values(model.indicator_constraint_info))
end

function MOI.get(
    model::Optimizer,
    ::MOI.ConstraintName,
    c::MOI.ConstraintIndex{<:MOI.VectorAffineFunction,<:MOI.Indicator},
)
    MOI.throw_if_not_valid(model, c)
    return _info(model, c).name
end

function MOI.set(
    model::Optimizer,
    ::MOI.ConstraintName,
    c::MOI.ConstraintIndex{<:MOI.VectorAffineFunction,S},
    name::String,
) where {S<:MOI.Indicator}
    MOI.throw_if_not_valid(model, c)
    _update_if_necessary(model)
    info = _info(model, c)
    info.name = name
    if !isempty(name)
        row = Cint(_info(model, c).row - 1)
        ret = GRBsetstrattrelement(model, "GenConstrName", row, name)
        _check_ret(model, ret)
    end
    return
end

function MOI.delete(
    model::Optimizer,
    c::MOI.ConstraintIndex{<:MOI.VectorAffineFunction,<:MOI.Indicator},
)
    MOI.throw_if_not_valid(model, c)
    row = _info(model, c).row
    ind = Ref{Cint}(row - 1)
    ret = GRBdelgenconstrs(model, 1, ind)
    _check_ret(model, ret)
    delete!(model.indicator_constraint_info, c.value)
    for info in values(model.indicator_constraint_info)
        if info.row > row
            info.row -= 1
        end
    end
    _require_update(model)
    return
end
