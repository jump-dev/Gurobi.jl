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
    return
end

end
