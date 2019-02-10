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

# Create the basic model. It is set up in such a way that the root relaxation
# is fractional.
@variable(model, 0 <= x <= 2, Int)
@variable(model, 0 <= y <= 4, Int)
@constraint(model, y <= 3.5 + x)
@constraint(model, y <= 4.1 - 0.2x)
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
    if cb_where == Gurobi.CB_MIPSOL
        # Gurobi has a feasible integer solution to the current problem.
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
    elseif cb_where == Gurobi.CB_MIPNODE
        # Gurobi has a fractional solution to the current problem.
        # Load the solution into the model.
        Gurobi.cbget_mipnode_rel(model, cb_data, cb_where)
        # Double check for sanity's sake that we have a feasible (given the
        # current constraints) point.
        @assert JuMP.primal_status(model) == MOI.FEASIBLE_POINT
        # Get the values of x and y.
        x_val, y_val = JuMP.value(x), JuMP.value(y)
        # Provide a heuristic solution. We don't need to provide a value for all
        # variables.
        Gurobi.cbsolution(model, cb_data, Dict(x => round(x_val)))
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
@test Gurobi.CB_POLLING in cb_calls    # A frequent callback to return control to the user
@test Gurobi.CB_PRESOLVE in cb_calls   # Entered presolve
@test !(Gurobi.CB_SIMPLEX in cb_calls) # Using the simplex method
@test Gurobi.CB_MIP in cb_calls        # Entered a MIP
@test Gurobi.CB_MIPSOL in cb_calls     # Found an integer solution
@test Gurobi.CB_MIPNODE in cb_calls    # Found a fractional solution
@test Gurobi.CB_MESSAGE in cb_calls    # Printing a message to the log
@test !(Gurobi.CB_BARRIER in cb_calls) # Using the barrier method
