# Gurobi.jl

[![Build Status](https://github.com/jump-dev/Gurobi.jl/workflows/CI/badge.svg?branch=master)](https://github.com/jump-dev/Gurobi.jl/actions?query=workflow%3ACI)
[![codecov](https://codecov.io/gh/jump-dev/Gurobi.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/jump-dev/Gurobi.jl)

[Gurobi.jl](https://github.com/jump-dev/Gurobi.jl) is a wrapper for the [Gurobi Optimizer](https://www.gurobi.com).

It has two components:
 - a thin wrapper around the complete C API
 - an interface to [MathOptInterface](https://github.com/jump-dev/MathOptInterface.jl)

The C API can be accessed via `Gurobi.GRBxx` functions, where the names and
arguments are identical to the C API. See the [Gurobi documentation](https://www.gurobi.com/documentation/9.0/refman/c_api_details.html)
for details.

## Affiliation

This wrapper is maintained by the JuMP community and is not officially
supported by Gurobi. However, we thank Gurobi for providing us with a license
to test Gurobi.jl on GitHub. If you are a commercial customer interested in
official support for Gurobi in Julia, let them know.

## License

`Gurobi.jl` is licensed under the [MIT License](https://github.com/jump-dev/Gurobi.jl/blob/master/LICENSE.md).

The underlying solver is a closed-source commercial product for which you must
[purchase a license](https://www.gurobi.com).

Free Gurobi licenses are available for [academics and students](https://www.gurobi.com/academia/academic-program-and-licenses/).

## Installation

First, obtain a license of Gurobi and install Gurobi solver, following the
instructions on [Gurobi's website](http://www.gurobi.com). Then, set the
`GUROBI_HOME` environment variable as appropriate and run `Pkg.add("Gurobi")`:

```julia
# On Windows, this might be
ENV["GUROBI_HOME"] = "C:\\Program Files\\gurobi1000\\win64"
# ... or perhaps ...
ENV["GUROBI_HOME"] = "C:\\gurobi1000\\win64"
# On Mac, this might be
ENV["GUROBI_HOME"] = "/Library/gurobi1000/mac64"

import Pkg
Pkg.add("Gurobi")
```
**Note: your path may differ. Check which folder you installed Gurobi in, and
update the path accordingly.**

By default, building Gurobi.jl will fail if the Gurobi library is not found.
This may not be desirable in certain cases, for example when part of a package's
test suite uses Gurobi as an optional test dependency, but Gurobi cannot be
installed on a CI server running the test suite. To support this use case, the
`GUROBI_JL_SKIP_LIB_CHECK` environment variable may be set (to any value) to
make Gurobi.jl installable (but not usable).

## Use with JuMP

To use Gurobi with [JuMP](https://github.com/jump-dev/JuMP.jl), use
`Gurobi.Optimizer`:
```julia
using JuMP, Gurobi
model = Model(Gurobi.Optimizer)
set_attribute(model, "TimeLimit", 100)
set_attribute(model, "Presolve", 0)
```

## MathOptInterface API

The Gurobi optimizer supports the following constraints and attributes.

List of supported objective functions:

 * [`MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}`](@ref)
 * [`MOI.ObjectiveFunction{MOI.ScalarQuadraticFunction{Float64}}`](@ref)
 * [`MOI.ObjectiveFunction{MOI.VariableIndex}`](@ref)
 * [`MOI.ObjectiveFunction{MOI.VectorAffineFunction{Float64}}`](@ref)

List of supported variable types:

 * [`MOI.Reals`](@ref)

List of supported constraint types:

 * [`MOI.ScalarAffineFunction{Float64}`](@ref) in [`MOI.EqualTo{Float64}`](@ref)
 * [`MOI.ScalarAffineFunction{Float64}`](@ref) in [`MOI.GreaterThan{Float64}`](@ref)
 * [`MOI.ScalarAffineFunction{Float64}`](@ref) in [`MOI.LessThan{Float64}`](@ref)
 * [`MOI.ScalarQuadraticFunction{Float64}`](@ref) in [`MOI.EqualTo{Float64}`](@ref)
 * [`MOI.ScalarQuadraticFunction{Float64}`](@ref) in [`MOI.GreaterThan{Float64}`](@ref)
 * [`MOI.ScalarQuadraticFunction{Float64}`](@ref) in [`MOI.LessThan{Float64}`](@ref)
 * [`MOI.VariableIndex`](@ref) in [`MOI.EqualTo{Float64}`](@ref)
 * [`MOI.VariableIndex`](@ref) in [`MOI.GreaterThan{Float64}`](@ref)
 * [`MOI.VariableIndex`](@ref) in [`MOI.Integer`](@ref)
 * [`MOI.VariableIndex`](@ref) in [`MOI.Interval{Float64}`](@ref)
 * [`MOI.VariableIndex`](@ref) in [`MOI.LessThan{Float64}`](@ref)
 * [`MOI.VariableIndex`](@ref) in [`MOI.Semicontinuous{Float64}`](@ref)
 * [`MOI.VariableIndex`](@ref) in [`MOI.Semiinteger{Float64}`](@ref)
 * [`MOI.VariableIndex`](@ref) in [`MOI.ZeroOne`](@ref)
 * [`MOI.VectorOfVariables`](@ref) in [`MOI.SOS1{Float64}`](@ref)
 * [`MOI.VectorOfVariables`](@ref) in [`MOI.SOS2{Float64}`](@ref)
 * [`MOI.VectorOfVariables`](@ref) in [`MOI.SecondOrderCone`](@ref)

List of supported model attributes:

 * [`MOI.HeuristicCallback()`](@ref)
 * [`MOI.LazyConstraintCallback()`](@ref)
 * [`MOI.Name()`](@ref)
 * [`MOI.ObjectiveSense()`](@ref)
 * [`MOI.UserCutCallback()`](@ref)

## Options

See the [Gurobi Documentation](https://www.gurobi.com/documentation/current/refman/parameters.html)
for a list and description of allowable parameters.

## Reusing the same Gurobi environment for multiple solves

When using this package via other packages such as [JuMP.jl](https://github.com/jump-dev/JuMP.jl),
the default behavior is to obtain a new Gurobi license token every time a model
is created. If you are using Gurobi in a setting where the number of concurrent
Gurobi uses is limited (for example, ["Single-Use" or "Floating-Use" licenses](http://www.gurobi.com/products/licensing-pricing/licensing-overview)),
you might instead prefer to obtain a single license token that is shared by all
models that your program solves. You can do this by passing a Gurobi Environment
object as the first parameter to `Gurobi.Optimizer`. For example, the follow
code snippet solves multiple problems with JuMP using the same license token:

```julia
using JuMP, Gurobi
const GRB_ENV = Gurobi.Env()

model1 = Model(() -> Gurobi.Optimizer(GRB_ENV))

# The solvers can have different options too
model2 = Model(() -> Gurobi.Optimizer(GRB_ENV))
set_attribute(model2, "OutputFlag", 0)
```

## Accessing Gurobi-specific attributes via JuMP

You can get and set Gurobi-specific variable, constraint, and model attributes
via JuMP as follows:

```julia
using JuMP, Gurobi
model = direct_model(Gurobi.Optimizer())
@variable(model, x >= 0)
@constraint(model, c, 2x >= 1)
@objective(model, Min, x)
MOI.set(model, Gurobi.ConstraintAttribute("Lazy"), c, 2)
optimize!(model)
MOI.get(model, Gurobi.VariableAttribute("LB"), x)  # Returns 0.0
MOI.get(model, Gurobi.ModelAttribute("NumConstrs")) # Returns 1
```

Note that we are using [JuMP in direct-mode](https://jump.dev/JuMP.jl/v0.20.0/solvers/#Direct-mode-1).

A complete list of supported Gurobi attributes can be found in
[their online documentation](https://www.gurobi.com/documentation/8.1/refman/attributes.html).

## Callbacks

Here is an example using Gurobi's solver-specific callbacks.

```julia
using JuMP, Gurobi, Test

model = direct_model(Gurobi.Optimizer())
@variable(model, 0 <= x <= 2.5, Int)
@variable(model, 0 <= y <= 2.5, Int)
@objective(model, Max, y)
cb_calls = Cint[]
function my_callback_function(cb_data, cb_where::Cint)
    # You can reference variables outside the function as normal
    push!(cb_calls, cb_where)
    # You can select where the callback is run
    if cb_where != GRB_CB_MIPSOL && cb_where != GRB_CB_MIPNODE
        return
    end
    # You can query a callback attribute using GRBcbget
    if cb_where == GRB_CB_MIPNODE
        resultP = Ref{Cint}()
        GRBcbget(cb_data, cb_where, GRB_CB_MIPNODE_STATUS, resultP)
        if resultP[] != GRB_OPTIMAL
            return  # Solution is something other than optimal.
        end
    end
    # Before querying `callback_value`, you must call:
    Gurobi.load_callback_variable_primal(cb_data, cb_where)
    x_val = callback_value(cb_data, x)
    y_val = callback_value(cb_data, y)
    # You can submit solver-independent MathOptInterface attributes such as
    # lazy constraints, user-cuts, and heuristic solutions.
    if y_val - x_val > 1 + 1e-6
        con = @build_constraint(y - x <= 1)
        MOI.submit(model, MOI.LazyConstraint(cb_data), con)
    elseif y_val + x_val > 3 + 1e-6
        con = @build_constraint(y + x <= 3)
        MOI.submit(model, MOI.LazyConstraint(cb_data), con)
    end
    if rand() < 0.1
        # You can terminate the callback as follows:
        GRBterminate(backend(model))
    end
    return
end
# You _must_ set this parameter if using lazy constraints.
MOI.set(model, MOI.RawOptimizerAttribute("LazyConstraints"), 1)
MOI.set(model, Gurobi.CallbackFunction(), my_callback_function)
optimize!(model)
@test termination_status(model) == MOI.OPTIMAL
@test primal_status(model) == MOI.FEASIBLE_POINT
@test value(x) == 1
@test value(y) == 2
```

See the [Gurobi documentation](https://www.gurobi.com/documentation/9.0/refman/cb_codes.html)
for other information that can be queried with `GRBcbget`.

### Common Performance Pitfall with JuMP

Gurobi's API works differently than most solvers. Any changes to the model are
not applied immediately, but instead go sit in a internal buffer (making any
modifications appear to be instantaneous) waiting for a call to [`GRBupdatemodel`](https://www.gurobi.com/documentation/9.0/refman/c_updatemodel.html)
(where the work is  done).

This leads to a common performance pitfall that has the following message as its
main symptom:

```
Warning: excessive time spent in model updates. Consider calling update less
frequently.
```

This often means the JuMP program was structured in such a way that Gurobi.jl
ends up calling `GRBupdatemodel` each iteration of a loop. Usually, it is
possible (and easy) to restructure the JuMP program in a way it stays
solver-agnostic and has a close-to-ideal performance with Gurobi.

To guide such restructuring it is good to keep in mind the following bits of
information:

1. `GRBupdatemodel` is only called if changes were done since last
   `GRBupdatemodel` (that is, if the internal buffer is not empty).
2. `GRBupdatemodel` is called when `JuMP.optimize!` is called, but this often is
not the source of the problem.
3. `GRBupdatemodel` *may* be called when *ANY* model attribute is queried *even
if that specific attribute was not changed*, and this often the source of the
problem.
4. The worst-case scenario is, therefore, a loop of modify-query-modify-query,
even if what is being modified and what is being queried are two completely
distinct things.

As an example, prefer:

```julia
# GOOD
model = Model(Gurobi.Optimizer)
@variable(model, x[1:100] >= 0)
# All modifications are done before any queries.
for i = 1:100
    set_upper_bound(x[i], i)
end
for i = 1:100
    # Only the first `lower_bound` query may trigger an `GRBupdatemodel`.
    println(lower_bound(x[i]))
end
```

to:

```julia
# BAD
model = Model(Gurobi.Optimizer)
@variable(model, x[1:100] >= 0)
for i = 1:100
    set_upper_bound(x[i], i)
    # `GRBupdatemodel` called on each iteration of this loop.
    println(lower_bound(x[i]))
end
```

## Common errors

### Using Gurobi v9.0 and you got an error like `Q not PSD`?

You need to set the NonConvex parameter:
```julia
model = Model(Gurobi.Optimizer)
set_optimizer_attribute(model, "NonConvex", 2)
```

### Gurobi Error 1009: Version number is XX.X, license is for version XX.X

First, please make sure that your license is correct for your Gurobi version.
See the [Gurobi documentation](https://support.gurobi.com/hc/en-us/articles/360034784572-How-do-I-check-for-a-valid-license-file-)
for details.

Once you are sure that the license and Gurobi versions match, re-install
Gurobi.jl by running:
```julia
import Pkg
Pkg.build("Gurobi")
```
