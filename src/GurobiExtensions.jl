module GurobiExtensions

using Gurobi, JuMP

export @lazy_constraint

macro lazy_constraint(model, cb_data, expr)
    code = quote
        lazy_con = @build_constraint $expr
        JuMP.backend($model)
        Gurobi.cblazy(
            JuMP.backend($model),
            $cb_data,
            JuMP.moi_function(lazy_con.func),
            lazy_con.set
        )
    end
    quote
        let
            $(esc(code))
        end
    end
end

function Gurobi.cbget_mipsol_sol(
        model::JuMP.Model, cb_data::Gurobi.CallbackData, cb_where::Int32)
    Gurobi.cbget_mipsol_sol(JuMP.backend(model), cb_data, cb_where)
end

function Gurobi.cbget_mipnode_rel(
        model::JuMP.Model, cb_data::Gurobi.CallbackData, cb_where::Int32)
    Gurobi.cbget_mipnode_rel(JuMP.backend(model), cb_data, cb_where)
end

function Gurobi.cbsolution(model::JuMP.Model, cb_data::Gurobi.CallbackData,
                           sol::Dict{JuMP.VariableRef, Float64})
    moi_sol = Dict{MOI.VariableIndex, Float64}(
        JuMP.index(variable) => value for (variable, value) in sol)
    Gurobi.cbsolution(JuMP.backend(model), cb_data, moi_sol)
end

end
