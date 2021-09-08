struct NumberOfObjectives <: MOI.AbstractModelAttribute end

function MOI.set(model::Optimizer, ::NumberOfObjectives, n::Integer)
    ret = GRBsetintattr(model, "NumObj", n)
    _check_ret(model, ret)
    return
end

function MOI.get(model::Optimizer, ::NumberOfObjectives)
    _update_if_necessary(model)
    n = Ref{Cint}()
    ret = GRBgetintattr(model, "NumObj", n)
    _check_ret(model, ret)
    return n[]
end

struct MultiObjectiveFunction <: MOI.AbstractModelAttribute
    index::Int
end

function MOI.set(
    model::Optimizer,
    attr::MultiObjectiveFunction,
    f::MOI.ScalarAffineFunction,
)
    num_vars = length(model.variable_info)
    obj = zeros(Float64, num_vars)
    for term in f.terms
        column = _info(model, term.variable).column
        obj[column] += term.coefficient
    end
    indices, coefficients = _indices_and_coefficients(model, f)
    _update_if_necessary(model)
    ret = GRBsetobjectiven(
        model,
        Cint(attr.index - 1),
        Cint(0),
        1.0,
        0.0,
        0.0,
        "Objective $(attr.index)",
        f.constant,
        length(indices),
        indices,
        coefficients,
    )
    _check_ret(model, ret)
    _require_update(model)
    return
end

struct MultiObjectivePriority <: MOI.AbstractModelAttribute
    index::Int
end

function MOI.set(
    model::Gurobi.Optimizer,
    attr::MultiObjectivePriority,
    priority::Int,
)
    env = GRBgetenv(model)
    ret = GRBsetintparam(env, "ObjNumber", attr.index - 1)
    _check_ret(env, ret)
    ret = GRBsetintattr(model, "ObjNPriority", priority)
    _check_ret(model, ret)
    _require_update(model)
    return
end

function MOI.get(model::Gurobi.Optimizer, attr::MultiObjectivePriority)
    _update_if_necessary(model)
    env = GRBgetenv(model)
    ret = GRBsetintparam(env, "ObjNumber", attr.index - 1)
    _check_ret(env, ret)
    n = Ref{Cint}()
    ret = GRBgetintattr(model, "ObjNPriority", n)
    _check_ret(model, ret)
    return Int(n[])
end

struct MultiObjectiveWeight <: MOI.AbstractModelAttribute
    index::Int
end

function MOI.set(
    model::Gurobi.Optimizer,
    attr::MultiObjectiveWeight,
    weight::Float64,
)
    env = GRBgetenv(model)
    ret = GRBsetintparam(env, "ObjNumber", attr.index - 1)
    _check_ret(env, ret)
    ret = GRBsetdblattr(model, "ObjNWeight", weight)
    _check_ret(model, ret)
    _require_update(model)
    return
end

function MOI.get(model::Gurobi.Optimizer, attr::MultiObjectiveWeight)
    _update_if_necessary(model)
    env = GRBgetenv(model)
    ret = GRBsetintparam(env, "ObjNumber", attr.index - 1)
    _check_ret(env, ret)
    weight = Ref{Cdouble}()
    ret = GRBgetdblattr(model, "ObjNWeight", weight)
    _check_ret(model, ret)
    return weight[]
end

struct MultiObjectiveValue <: MOI.AbstractModelAttribute
    index::Int
end

function MOI.get(model::Gurobi.Optimizer, attr::MultiObjectiveValue)
    env = GRBgetenv(model)
    ret = GRBsetintparam(env, "ObjNumber", attr.index - 1)
    _check_ret(env, ret)
    val = Ref{Cdouble}()
    ret = GRBgetdblattr(model, "ObjNVal", val)
    _check_ret(model, ret)
    return val[]
end
