# Copyright (c) 2015 Dahua Lin, Miles Lubin, Joey Huchette, Iain Dunning, and
#   contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

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

"""
    MultiObjectiveAttribute(index::Int, name::String)

A Gurobi-specific attribute for interacting with multi-objective attributes.

`index` is the 1-based index of the objective to query. It corresponds to
setting `ObjNumber` to `index - 1`.

See also: `MultiObjectiveWeight`, `MultiObjectivePriority`, `MultiObjectiveValue`.

## Examples

```julia
MOI.get(model, MultiObjectiveAttribute(1, "ObjNPriority"))
MOI.get(model, MultiObjectiveAttribute(2, "ObjNVal"))
MOI.set(model, MultiObjectiveAttribute(2, "ObjNWeight"), 2.0)
```
"""
struct MultiObjectiveAttribute <: MOI.AbstractModelAttribute
    index::Int
    name::String
end

function MOI.set(model::Optimizer, attr::MultiObjectiveAttribute, value)
    env = GRBgetenv(model)
    ret = GRBsetintparam(env, "ObjNumber", attr.index - 1)
    _check_ret(env, ret)
    MOI.set(model, ModelAttribute(attr.name), value)
    return
end

function MOI.get(model::Optimizer, attr::MultiObjectiveAttribute)
    _update_if_necessary(model)
    env = GRBgetenv(model)
    ret = GRBsetintparam(env, "ObjNumber", attr.index - 1)
    _check_ret(env, ret)
    return MOI.get(model, ModelAttribute(attr.name))
end

struct MultiObjectivePriority <: MOI.AbstractModelAttribute
    index::Int
end

function MOI.set(model::Optimizer, attr::MultiObjectivePriority, value::Int)
    MOI.set(model, MultiObjectiveAttribute(attr.index, "ObjNPriority"), value)
    return
end

function MOI.get(model::Optimizer, attr::MultiObjectivePriority)
    return MOI.get(model, MultiObjectiveAttribute(attr.index, "ObjNPriority"))
end

struct MultiObjectiveWeight <: MOI.AbstractModelAttribute
    index::Int
end

function MOI.set(model::Optimizer, attr::MultiObjectiveWeight, weight::Float64)
    MOI.set(model, MultiObjectiveAttribute(attr.index, "ObjNWeight"), weight)
    return
end

function MOI.get(model::Optimizer, attr::MultiObjectiveWeight)
    return MOI.get(model, MultiObjectiveAttribute(attr.index, "ObjNWeight"))
end

struct MultiObjectiveValue <: MOI.AbstractModelAttribute
    index::Int
end

function MOI.get(model::Optimizer, attr::MultiObjectiveValue)
    return MOI.get(model, MultiObjectiveAttribute(attr.index, "ObjNVal"))
end

function MOI.set(
    model::Optimizer,
    ::MOI.ObjectiveFunction{F},
    f::F,
) where {F<:MOI.VectorAffineFunction{Float64}}
    MOI.set(model, NumberOfObjectives(), MOI.output_dimension(f))
    for (i, fi) in enumerate(MOI.Utilities.eachscalar(f))
        MOI.set(model, MultiObjectiveFunction(i), fi)
    end
    model.objective_type = _VECTOR_AFFINE
    return
end

function MOI.get(
    model::Optimizer,
    ::MOI.ObjectiveFunction{MOI.VectorAffineFunction{Float64}},
)
    env = GRBgetenv(model)
    F = MOI.ScalarAffineFunction{Float64}
    f = F[]
    for i in 1:MOI.get(model, NumberOfObjectives())
        ret = GRBsetintparam(env, "ObjNumber", i - 1)
        _check_ret(env, ret)
        push!(f, _get_affine_objective(model; is_multiobjective = true))
    end
    return MOI.Utilities.operate(vcat, Float64, f...)
end

function MOI.supports(
    model::Optimizer,
    ::MOI.ObjectiveFunction{MOI.VectorAffineFunction{Float64}},
)
    return true
end
