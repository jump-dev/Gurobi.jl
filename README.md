**Gurobi.jl underwent a major rewrite between versions 0.8.1 and 0.9.0. Users of
JuMP should see no breaking changes, but if you used the lower-level C API
(e.g., for callbacks), you will need to update your code accordingly. For a full
description of the changes, read [this discourse post](https://discourse.julialang.org/t/ann-upcoming-breaking-changes-to-cplex-jl-and-gurobi-jl/47814).**

**To revert to the old API, use:**
```julia
import Pkg
Pkg.add(Pkg.PackageSpec(name = "Gurobi", version = v"0.8"))
```
**Then restart Julia for the change to take effect.**

# Gurobi.jl

[![Build Status](https://github.com/jump-dev/Gurobi.jl/workflows/CI/badge.svg?branch=master)](https://github.com/jump-dev/Gurobi.jl/actions?query=workflow%3ACI)
[![codecov](https://codecov.io/gh/jump-dev/Gurobi.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/jump-dev/Gurobi.jl)

Gurobi.jl is a wrapper for the [Gurobi Optimizer](https://www.gurobi.com).

It has two components:
 - a thin wrapper around the complete C API
 - an interface to [MathOptInterface](https://github.com/jump-dev/MathOptInterface.jl)

The C API can be accessed via `Gurobi.GRBxx` functions, where the names and
arguments are identical to the C API. See the [Gurobi documentation](https://www.gurobi.com/documentation/9.0/refman/c_api_details.html)
for details.

*Note: This wrapper is maintained by the JuMP community and is not officially
supported by Gurobi. However, we thank Gurobi for providing us with a license
to test Gurobi.jl on GitHub. If you are a commercial customer interested in
official support for Gurobi in Julia, let them know!.*

## Installation

**Minimum version requirement:** Gurobi.jl requires Gurobi version 9.0 or 9.1 or 9.5.

First, obtain a license of Gurobi and install Gurobi solver, following the
instructions on [Gurobi's website](http://www.gurobi.com). Then, set the
`GUROBI_HOME` environment variable as appropriate and run `Pkg.add("Gurobi")`,
the `Pkg.build("Gurobi")`. For example:
```julia
# On Windows, this might be
ENV["GUROBI_HOME"] = "C:\\Program Files\\gurobi950\\win64"
# ... or perhaps ...
ENV["GUROBI_HOME"] = "C:\\gurobi950\\win64"
import Pkg
Pkg.add("Gurobi")
Pkg.build("Gurobi")

# On Mac, this might be
ENV["GUROBI_HOME"] = "/Library/gurobi950/mac64"
import Pkg
Pkg.add("Gurobi")
Pkg.build("Gurobi")
```
**Note: your path may differ. Check which folder you installed Gurobi in, and
update the path accordingly.**

By default, `build`ing *Gurobi.jl* will fail if the Gurobi library is not found.
This may not be desirable in certain cases, for example when part of a package's
test suite uses Gurobi as an optional test dependency, but Gurobi cannot be
installed on a CI server running the test suite. To support this use case, the
`GUROBI_JL_SKIP_LIB_CHECK` environment variable may be set (to any value) to
make *Gurobi.jl* installable (but not usable).

## Use with JuMP

We highly recommend that you use the *Gurobi.jl* package with higher level
packages such as [JuMP.jl](https://github.com/jump-dev/JuMP.jl).

This can be done using the ``Gurobi.Optimizer`` object. Here is how to create a
*JuMP* model that uses Gurobi as the solver.
```julia
using JuMP, Gurobi

model = Model(Gurobi.Optimizer)
set_optimizer_attribute(model, "TimeLimit", 100)
set_optimizer_attribute(model, "Presolve", 0)
```
See the [Gurobi Documentation](https://www.gurobi.com/documentation/current/refman/parameters.html)
for a list and description of allowable parameters.

## Reusing the same Gurobi environment for multiple solves

When using this package via other packages such as [JuMP.jl](https://github.com/jump-dev/JuMP.jl),
the default behavior is to obtain a new Gurobi license token every time a model
is created. If you are using Gurobi in a setting where the number of concurrent
Gurobi uses is limited (e.g. ["Single-Use" or "Floating-Use" licenses](http://www.gurobi.com/products/licensing-pricing/licensing-overview)),
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
set_optimizer_attribute(model2, "OutputFlag", 0)
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
MOI.set(model, MOI.RawParameter("LazyConstraints"), 1)
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
`GRBupdatemodel` (i.e., the internal buffer is not empty).
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

## Using Gurobi v9.0 and you got an error like `Q not PSD`?

You need to set the NonConvex parameter:
```julia
model = Model(Gurobi.Optimizer)
set_optimizer_attribute(model, "NonConvex", 2)
```
