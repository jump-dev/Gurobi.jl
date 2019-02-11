module GurobiJuMPExtensions

using Gurobi, JuMP

function JuMP._set_callbacks(optimizer::Gurobi.Optimizer,
                             lazy_callback, heuristic_callback)
    MOI.set(optimizer, Gurobi.CallbackFunction(), (cb_data, cb_where) -> begin
        if lazy_callback !== nothing && cb_where == Gurobi.CB_MIPSOL
            Gurobi.cbget_mipsol_sol(optimizer, cb_data, cb_where)
            lazy_callback(cb_data)
        elseif heuristic_callback !== nothing && cb_where == Gurobi.CB_MIPNODE
            Gurobi.cbget_mipnode_rel(optimizer, cb_data, cb_where)
            heuristic_callback(cb_data)
        end
        return
    end)
end

function JuMP.add_lazy_constraint(
        model::Gurobi.Optimizer, cb_data::Gurobi.CallbackData, func, set)
    Gurobi.cblazy(model, cb_data, func, set)
end

function JuMP.add_heuristic_solution(model::Gurobi.Optimizer,
                                     cb_data::Gurobi.CallbackData,
                                     sol::Dict{JuMP.VariableRef, Float64})
    Gurobi.cbsolution(model, cb_data, Dict(
        JuMP.index(variable) => value for (variable, value) in sol))
end

end
