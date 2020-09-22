# Gurobi.jl

The [Gurobi](http://www.gurobi.com) Optimizer is a commercial optimization
solver for a variety of mathematical programming problems, including linear
programming (LP), quadratic programming (QP), quadratically constrained
programming (QCP), mixed integer linear programming (MILP), mixed-integer
quadratic programming (MIQP), and mixed-integer quadratically constrained
programming (MIQCP).

*The Gurobi wrapper for Julia is community driven and not officially supported
by Gurobi. If you are a commercial customer interested in official support for
Julia from Gurobi, let them know!*

## Installation

First, obtain a license of Gurobi and install Gurobi solver, following the
instructions on [Gurobi's website](http://www.gurobi.com).

   **The minimum version supported by *Gurobi.jl* is Gurobi v9.0.**

Then, run the following:

```julia
import Pkg
# On Mac
ENV["GUROBI_HOME"] = "/Library/gurobi902/mac64"
# On Windows
ENV["GUROBI_HOME"] = "C:/Program Files/gurobi902/win64"
# ... or perhaps ...
ENV["GUROBI_HOME"] = "C:/gurobi902/win64"

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
