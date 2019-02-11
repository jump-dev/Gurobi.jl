module GurobiJuMPExtensions

using Gurobi, JuMP

function JuMP._set_callbacks(model::Model, optimizer::Gurobi.Optimizer,
                             lazy_callback, heuristic_callback)
    MOI.set(optimizer, Gurobi.CallbackFunction(), (cb_data, cb_where) -> begin
        if lazy_callback !== nothing && cb_where == Gurobi.CB_MIPSOL
            Gurobi.cbget_mipsol_sol(optimizer, cb_data, cb_where)
            lazy_callback(model, cb_data)
        elseif heuristic_callback !== nothing && cb_where == Gurobi.CB_MIPNODE
            Gurobi.cbget_mipnode_rel(optimizer, cb_data, cb_where)
            heuristic_callback(model, cb_data)
        end
        return
    end)
end

function JuMP.add_lazy_constraint(
        model::Model, cb_data::Gurobi.CallbackData, func, set)
    Gurobi.cblazy(JuMP.backend(model), cb_data, func, set)
end

function JuMP.add_heuristic_solution(model::Model, cb_data::Gurobi.CallbackData,
                                     sol::Dict{JuMP.VariableRef, Float64})
    Gurobi.cbsolution(JuMP.backend(model), cb_data, Dict(
        JuMP.index(variable) => value for (variable, value) in sol))
end

end
