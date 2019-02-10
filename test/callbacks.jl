using JuMP, Gurobi, Test

# This special macro loads a module called GurobiExtensions, which exports the
# `@lazy_constraint` macro.
Gurobi.@load_extensions

# Create a JuMP model using Gurobi in direct-mode. We have to choose some
# settings like `Presolve=0` to prevent Gurobi from trivially solving this
# problem. Since we're using lazy constraints, we also have to set
# `LazyConstraints=1`.
model = JuMP.direct_model(
    Gurobi.Optimizer(
        OutputFlag=1, Cuts=0, Presolve=0, Heuristics=0, LazyConstraints=1
    )
)

# Now add the decision variables and the objective.
@variable(model, 0 <= x <= 2, Int)
@variable(model, 0 <= y <= 2, Int)
@objective(model, Max, y)

# This vector is going to cache the value of `cb_where` everytime our callback
# gets called.
cb_calls = Int32[]

# Here is the callback function. We set the solver-specific attribute
# `Gurobi.CallbackFunction()`, passing a function (in this case, anonymous) as
# the third argument.
MOI.set(model, Gurobi.CallbackFunction(), (cb_data, cb_where) -> begin
    # Cache the value of `cb_where` in `cb_calls`. Note how we can access
    # variables in the outer scope.
    push!(cb_calls, cb_where)
    # Check where this callback is being called from. We can only add lazy
    # constraints when `cb_where == CB_MIPSOL`.
    if cb_where == Gurobi.CB_MIPSOL
        # Load the intger solution into the model. This lets us query
        # `JuMP.value` on variables.
        Gurobi.cbget_mipsol_sol(model, cb_data, cb_where)
        # Double check for sanity's sake that we have a feasible (given the
        # current constraints) point.
        @assert JuMP.primal_status(model) == MOI.FEASIBLE_POINT
        # Get the values of x and y.
        x_val, y_val = JuMP.value(x), JuMP.value(y)
        # Add the lazy constraints using the `@lazy_constraint` macro that was
        # loaded by `Gurobi.@load_extensions`.
        if y_val - x_val > 1 + 1e-6
            @lazy_constraint(model, cb_data, y <= 1 + x)
        elseif y_val + x_val > 3 + 1e-6
            @lazy_constraint(model, cb_data, y <= 3 - x)
        end
    end
    return
end)

# Solve the model.
JuMP.optimize!(model)

# Check the solution.
@test JuMP.value(x) == 1
@test JuMP.value(y) == 2

# Check that our callback has been used.
@test length(cb_calls) > 0
@test Gurobi.CB_MESSAGE in cb_calls
@test Gurobi.CB_PRESOLVE in cb_calls
@test Gurobi.CB_MIPSOL in cb_calls
