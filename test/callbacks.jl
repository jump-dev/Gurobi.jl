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
cb_calls = Symbol[]

function my_lazy_callback(cb_data)
    push!(cb_calls, :lazy)
    # Double check for sanity's sake that we have a feasible (given the
    # current constraints) point.
    @assert primal_status(model) == MOI.FEASIBLE_POINT
    # Get the values of x and y.
    x_val, y_val = value(x), value(y)
    # Add the lazy constraints using the `@lazy_constraint` macro that was
    # loaded by `@load_extensions`.
    if y_val - x_val > 1.1 + 1e-6
        @lazy_constraint(model, cb_data, y <= 1.1 + x)
    elseif y_val + x_val > 3 + 1e-6
        @lazy_constraint(model, cb_data, y <= 3 - x)
    end
end

function my_heuristic_callback(cb_data)
    push!(cb_calls, :heuristic)
    # Double check for sanity's sake that we have a feasible (given the
    # current constraints) point.
    @assert primal_status(model) == MOI.FEASIBLE_POINT
    # Get the values of x and y.
    x_val, y_val = value(x), value(y)
    # Provide a heuristic solution.
    add_heuristic_solution(
        model, cb_data, Dict(x => 1.0, y => 2.0))
end

set_callbacks(model;
    lazy = my_lazy_callback,
    heuristic = my_heuristic_callback
)

# Solve the model.
optimize!(model)

# Check the solution.
@test value(x) == 1
@test value(y) == 2

# Check that our callback has been used.
@test length(cb_calls) > 0
@test :lazy in cb_calls
@test :heuristic in cb_calls
