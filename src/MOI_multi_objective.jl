struct NumberOfObjectives <: MOI.AbstractModelAttribute end

function MOI.set(model::Optimizer, ::NumberOfObjectives, n::Integer)
    set_intattr!(model.inner, "NumObj", n)
    return
end

function MOI.get(model::Optimizer, ::NumberOfObjectives)
    return get_intattr(model.inner, "NumObj")
end

struct MultiObjectiveFunction <: MOI.AbstractModelAttribute
    index::Int
end

function c_api_setobjectiven(
    model,
    index::Cint,
    priority::Cint,
    weight::Cdouble,
    abstol::Cdouble,
    reltol::Cdouble,
    name::String,
    constant::Cdouble,
    lnz::Cint,
    lind::Vector{Cint},
    lval::Vector{Cdouble}
)
    ret = @grb_ccall(
        setobjectiven,
        Cint,
        (
            Ptr{Cvoid}, Cint, Cint, Cdouble, Cdouble, Cdouble, Ptr{Cchar},
            Cdouble, Cint, Ptr{Cint}, Ptr{Cdouble}
        ),
        model, index, priority, weight, abstol, reltol, name, constant, lnz, lind, lval
    )
    if ret != 0
        throw(GurobiError(model.env, ret))
    end
end

function MOI.set(
    model::Optimizer,
    attr::MultiObjectiveFunction,
    f::MOI.ScalarAffineFunction
)
    num_vars = length(model.variable_info)
    obj = zeros(Float64, num_vars)
    for term in f.terms
        column = _info(model, term.variable_index).column
        obj[column] += term.coefficient
    end
    indices, coefficients = _indices_and_coefficients(model, f)
    _update_if_necessary(model)
    c_api_setobjectiven(
        model.inner,
        Cint(attr.index - 1),
        Cint(0),
        1.0,
        0.0,
        0.0,
        "",
        f.constant,
        Cint(length(indices)),
        Cint.(indices) .- Cint(1),
        coefficients
    )
    _require_update(model)
    return
end

struct MultiObjectivePriority <: MOI.AbstractModelAttribute
    index::Int
end

function MOI.set(
    model::Gurobi.Optimizer, attr::MultiObjectivePriority, priority::Int
)
    set_int_param!(model.inner, "ObjNumber", attr.index - 1)
    set_intattr!(model.inner, "ObjNPriority", priority)
    _require_update(model)
    return
end

function MOI.get(
    model::Gurobi.Optimizer, attr::MultiObjectivePriority
)
    _update_if_necessary(model)
    set_int_param!(model.inner, "ObjNumber", attr.index - 1)
    return get_intattr(model.inner, "ObjNPriority")
end

struct MultiObjectiveWeight <: MOI.AbstractModelAttribute
    index::Int
end

function MOI.set(
    model::Gurobi.Optimizer, attr::MultiObjectiveWeight, weight::Float64
)
    set_int_param!(model.inner, "ObjNumber", attr.index - 1)
    set_dblattr!(model.inner, "ObjNWeight", weight)
    _require_update(model)
    return
end

function MOI.get(
    model::Gurobi.Optimizer, attr::MultiObjectiveWeight
)
    _update_if_necessary(model)
    set_int_param!(model.inner, "ObjNumber", attr.index - 1)
    return get_dblattr(model.inner, "ObjNWeight")
end

struct MultiObjectiveValue <: MOI.AbstractModelAttribute
    index::Int
end

function MOI.get(
    model::Gurobi.Optimizer, attr::MultiObjectiveValue
)
    set_int_param!(model.inner, "ObjNumber", attr.index - 1)
    return get_dblattr(model.inner, "ObjNVal")
end
