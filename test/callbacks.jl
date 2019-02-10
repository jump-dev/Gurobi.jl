using JuMP, Gurobi, Test

model = JuMP.direct_model(Gurobi.Optimizer(
    OutputFlag=1, Cuts=0, Presolve=0, Heuristics=0, LazyConstraints=1))
@variable(model, 0 <= x <= 2, Int)
@variable(model, 0 <= y <= 2, Int)
@objective(model, Max, y)

cb_calls = Int32[]
function callback_function(cb_data, cb_where)
    push!(cb_calls, cb_where)
    if cb_where == Gurobi.CB_MIPSOL
        Gurobi.cbget_mipsol_sol(JuMP.backend(model), cb_data, cb_where)
        x_val, y_val = JuMP.value(x), JuMP.value(y)
        if y_val - x_val > 1 + 1e-6
            @lazy_constraint(model, cb_data, y <= 1 + x)
        elseif y_val + x_val > 3 + 1e-6
            @lazy_constraint(model, cb_data, y <= 3 - x)
        end
    end
    return
end
MOI.set(model, Gurobi.CallbackFunction(), callback_function)

JuMP.optimize!(model)

@test JuMP.value(x) == 1
@test JuMP.value(y) == 2
@test length(cb_calls) > 0
@test Gurobi.CB_MESSAGE in cb_calls
@test Gurobi.CB_PRESOLVE in cb_calls
@test Gurobi.CB_MIPSOL in cb_calls

macro lazy_constraint(model, cb_data, expr)
    quote
        lazy_con = JuMP.@build_constraint($(expr))
        Gurobi.cblazy(
            JuMP.backend($(esc(model))), $(esc(cb_data)),
            JuMP.moi_function(lazy_con.func), lazy_con.set)
    end
end
