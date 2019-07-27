# Gurobi.jl

The [Gurobi](http://www.gurobi.com) Optimizer is a commercial optimization solver for a variety of mathematical programming problems, including linear programming (LP), quadratic programming (QP), quadratically constrained programming (QCP), mixed integer linear programming (MILP), mixed-integer quadratic programming (MIQP), and mixed-integer quadratically constrained programming (MIQCP).

*The Gurobi wrapper for Julia is community driven and not officially supported by Gurobi. If you are a commercial customer interested in official support for Julia from Gurobi, let them know!*

## Use with JuMP

We highly recommend that you use the *Gurobi.jl* package with higher level packages such as [JuMP.jl](https://github.com/JuliaOpt/JuMP.jl).

This can be done using the ``Gurobi.Optimizer`` object. Here is how to create a *JuMP* model that uses Gurobi as the solver. Parameters are passed as keyword arguments:
```julia
using JuMP, Gurobi

model = Model(with_optimizer(Gurobi.Optimizer, Presolve=0, OutputFlag=0))
```
See the [Gurobi Documentation](https://www.gurobi.com/documentation/8.1/refman/parameters.html) for a list and description of allowable parameters.

## Installation

Here is the procedure to setup this package:

1. Obtain a license of Gurobi and install Gurobi solver, following the instructions on [Gurobi's website](http://www.gurobi.com).

   **The minimum version supported by *Gurobi.jl* is Gurobi v7.0.**

2. Make sure the ``GUROBI_HOME`` environmental variable is set to the path of the Gurobi directory. This is part of a standard installation. The Gurobi library will be searched for in ``GUROBI_HOME/lib`` on unix platforms and ``GUROBI_HOME\bin`` on Windows. If the library is not found, check that your version is listed in ``deps/build.jl``.

3. Install this package using ``Pkg.add("Gurobi")``.

4. Now, you can start using it.

By default, `build`ing *Gurobi.jl* will fail if the Gurobi library is not found. This may not be desirable in certain cases, for example when part of a package's test suite uses Gurobi as an optional test dependency, but Gurobi cannot be installed on a CI server running the test suite. To support this use case, the `GUROBI_JL_SKIP_LIB_CHECK` environment variable may be set (to any value) to make *Gurobi.jl* installable (but not usable).

## "LoadError: Unable to locate Gurobi installation"

- Make sure that you have downloaded and installed Gurobi from [gurobi.com](https://gurobi.com). Also, make sure that you validate the license.

- Check that you have a Gurobi version between 7.0 and 8.1.

- Make sure that the `GUROBI_HOME` environment variable is set correctly. You can see the current value as follows
```julia
julia> ENV["GUROBI_HOME"]
"C:\\gurobi801\\win64"
```
If it is not set correctly (e.g., you get an error `Key "GUROBI_HOME" not found`), you can set it as follows
```julia
julia> ENV["GUROBI_HOME"] = "/replace/this/with/the/path/to/gurobi"

julia> import Pkg; Pkg.build("Gurobi")
```
The Gurobi library (`gurobiXXX.dll` on Windows, `gurobiXXX.so` on Unix, and `gurobiXXX.dylib` in OSX where `XXX` is a version) will be searched for in ``GUROBI_HOME/lib`` on unix platforms and ``GUROBI_HOME\bin`` on Windows.

*Most users should not need to use the low-level API detailed in the following sections.*

## API Overview

This package provides both APIs at different levels for constructing models and solving optimization problems.

#### Gurobi Environment

A Gurobi model is always associated with an Gurobi environment, which maintains a solver configuration. By setting parameters to this environment, one can control or tune the behavior of a Gurobi solver.

To construct a Gurobi Environment, one can write:
```
env = Gurobi.Env()
```

This package provides functions to get and set parameters:

```julia
getparam(env, name)       # get the value of a parameter
setparam!(env, name, v)   # set the value of a parameter
setparams!(env, name1=value1, name2=value2, ...)  # set parameters using keyword arguments
```

You may refer to Gurobi's [Parameter Reference](http://www.gurobi.com/documentation/8.1/refman/parameters.html) for the whole list of parameters.

Here are some simple examples
```julia
setparam!(env, "Method", 2)   # choose to use Barrier method
setparams!(env; IterationLimit=100, Method=1) # set the maximum iterations and choose to use Simplex method
```

These parameters may be used directly with the ``GurobiSolver`` object used by MathProgBase. For example:
```julia
solver = GurobiSolver(Method=2)
solver = GurobiSolver(Method=1, IterationLimit=100.)
```

#### High-level API

If the objective coefficients and the constraints have already been given, one may use a high-level function ``gurobi_model`` to construct a model:

```julia
gurobi_model(env, ...)
```
One can use keyword arguments to specify the models:
* ``name``:  the model name.
* ``sense``: the sense of optimization (a symbol, which can be either ``:minimize`` (default) or ``:maximize``).
* ``f``:   the linear coefficient vector.
* ``H``:   the quadratic coefficient matrix (can be dense or sparse).
* ``A``:   the coefficient matrix of the linear inequality constraints.
* ``b``:   the right-hand-side of the linear inequality constraints.
* ``Aeq``:  the coefficient matrix of the equality constraints.
* ``beq``:  the right-hand-side of the equality constraints.
* ``lb``:   the variable lower bounds.
* ``ub``:   the variable upper bounds.

This function constructs a model that represents the following problem:
```
objective:  (1/2) x' H x + f' x

      s.t.   A x <= b
           Aeq x = beq
         lb <= x <= ub
```

The caller *must* specify ``f`` using a non-empty vector, while other keyword arguments are optional. When ``H`` is omitted, this reduces to an LP problem. When ``lb`` is omitted, the variables are not lower bounded, and when ``ub`` is omitted, the variables are not upper bounded.


#### Low-level API

This package also provides functions to build the model from scratch and gradually add variables and constraints.
To construct an empty model, one can write:
```julia
env = Gurobi.Env()    # creates a Gurobi environment

model = Gurobi.Model(env, name)   # creates an empty model
model = Gurobi.Model(env, name, sense)
```

Here, ``sense`` is a symbol, which can be either ``:minimize`` or ``:maximize`` (default to ``:minimize`` when omitted).

Then, the following functions can be used to add variables and constraints to the model:
```julia
## add variables

add_var!(model, vtype, c)   # add an variable with coefficient c
                            # vtype can be either of
                            # - GRB_CONTINUOUS  (for continuous variable)
                            # - GRB_INTEGER (for integer variable)
                            # - GRB_BINARY (for binary variable, i.e. 0/1)

add_cvar!(model, c)            # add a continuous variable
add_cvar!(model, c, lb, ub)    # add a continuous variable with specified bounds

add_ivar!(model, c)            # add an integer variable
add_ivar!(model, c, lb, ub)    # add an integer variable with specified bounds

add_bvar!(model, c)            # add a binary variable

## add constraints

# add a constraint with non-zero coefficients on specific variables.
# rel can be '<', '>', or '='
add_constr!(model, inds, coeffs, rel, rhs)

# add a constraint with coefficient vector for all variables.
add_constr!(model, coeffs, rel, rhs)

# add constraints using CSR format
add_constrs!(model, cbegin, inds, coeffs, rel, rhs)

# add constraints using a matrix: A x (rel) rhs
add_constrs!(model, A, rel, rhs)  # here A can be dense or sparse

# add constraints using a transposed matrix: At' x (rel) rhs
# this is usually more efficient than add_constrs!
add_constrs_t!(model, At, rel, rhs)  # here At can be dense or sparse

# add a range constraint
add_rangeconstr!(model, inds, coeffs, lb, ub)

# add range constraints using CSR format
add_rangeconstrs!(model, cbegin, inds, coeffs, lb, ub)

# add range constraints using a matrix:  lb <= A x <= ub
add_rangeconstrs!(model, A, lb, ub)  # here A can be dense or sparse

# add range constraints using a transposed matrix: lb <= At' x <= ub
# this is usually more efficient than add_rangeconstrs!
add_rangeconstrs_t!(model, At, lb, ub)  # here At can be dense or sparse
```

#### Modify Problem

It is not uncommon in practice that one would like to adjust the objective coefficients and solve the problem again. This package provides a function ``set_objcoeffs!`` for this purpose:
```julia
set_objcoeffs!(model, new_coeffs)
 # ... one can also call add_constr! and friends to add additional constraints ...
update_model!(model)   # changes take effect after this
optimize(model)
```


## Examples

The usage of this package is straight forward. Below, we use several examples to demonstrate the use of this package to solve optimization problems.

### Linear Programming Examples

Problem formulation:
```
maximize x + y

s.t. 50 x + 24 y <= 2400
     30 x + 33 y <= 2100
     x >= 45, y >= 5
```

Below, we show how this problem can be constructed and solved in different ways.

##### Example 1.1: High-level Linear Programming API

Using the ``gurobi_model`` function:

```julia
using Gurobi

env = Gurobi.Env()

# set presolve to 0
setparam!(env, "Presolve", 0)

 # construct the model
model = gurobi_model(env;
    name = "lp_01",
    f = ones(2),
    A = [50. 24.; 30. 33.],
    b = [2400., 2100.],
    lb = [5., 45.])

 # run optimization
optimize(model)

 # show results
sol = get_solution(model)
println("soln = $(sol)")

objv = get_objval(model)
println("objv = $(objv)")
```

##### Example 1.2: Low-level Linear Programming API

```julia
using Gurobi

env = Gurobi.Env()

# set presolve to 0
setparam!(env, "Presolve", 0)

 # creates an empty model ("lp_01" is the model name)
model = Gurobi.Model(env, "lp_01", :maximize)

 # add variables
 # add_cvar!(model, obj_coef, lower_bound, upper_bound)
add_cvar!(model, 1.0, 45., Inf)  # x: x >= 45
add_cvar!(model, 1.0,  5., Inf)  # y: y >= 5

 # For Gurobi, you have to call update_model to have the
 # lastest changes take effect
update_model!(model)

 # add constraints
 # add_constr!(model, coefs, sense, rhs)
add_constr!(model, [50., 24.], '<', 2400.) # 50 x + 24 y <= 2400
add_constr!(model, [30., 33.], '<', 2100.) # 30 x + 33 y <= 2100
update_model!(model)

println(model)

 # perform optimization
optimize(model)
```

You may also add variables and constraints in batch, as:

```julia
 # add mutliple variables in batch
add_cvars!(model, [1., 1.], [45., 5.], Inf)

 # add multiple constraints in batch
A = [50. 24.; 30. 33.]
b = [2400., 2100.]
add_constrs!(model, A, '<', b)
```

##### Example 1.3: Linear programming (MATLAB-like style)

You may also specify and solve the entire problem in one function call, using the
solver-independent [MathProgBase](https://github.com/JuliaOpt/MathProgBase.jl) package.

Julia code:
```julia
using MathProgBase, Gurobi

f = [1., 1.]
A = [50. 24.; 30. 33.]
b = [2400., 2100.]
lb = [5., 45.]

# pass params as keyword arguments to GurobiSolver
solution = linprog(f, A, '<', b, lb, Inf, GurobiSolver(Presolve=0))
```

##### Example 1.4: Linear programming with JuMP (Algebraic model)

Using [JuMP](https://github.com/JuliaOpt/JuMP.jl), we can specify linear programming problems using a more
natural algebraic approach.

```julia
using JuMP, Gurobi

# pass params as keyword arguments to GurobiSolver
model = Model(with_optimizer(Gurobi.Optimizer, Presolve=0))

@variable(model, x >= 5)
@variable(model, y >= 45)

@objective(model, Min, x + y)
@constraint(model, 50x + 24y <= 2400)
@constraint(model, 30x + 33y <= 2100)

optimize!(model)
println("Optimal objective: ", objective_value(model),
	". x = ", value(x), " y = ", value(y))
```

### Quadratic programming Examples

Problem formulation:
```
minimize x^2 + xy + y^2 + yz + z^2

s.t.  x + 2 y + 3 z >= 4
      x +   y       >= 1
```

##### Example 2.1: High-level Quadratic Programming API

using the function ``gurobi_model``:
```julia
using Gurobi

env = Gurobi.Env()

model = gurobi_model(env;
        name = "qp_01",
        H = [2. 1. 0.; 1. 2. 1.; 0. 1. 2.],
        f = [0., 0., 0.],
        A = -[1. 2. 3.; 1. 1. 0.],
        b = -[4., 1.])
optimize(model)
```

##### Example 2.2: Low-level Quadratic Programming API

```julia
using Gurobi

env = Gurobi.Env()

model = Gurobi.Model(env, "qp_01")

add_cvars!(model, [1., 1.], 0., Inf)
update_model!(model)

 # add quadratic terms: x^2, x * y, y^2
 # add_qpterms!(model, rowinds, colinds, coeffs)
add_qpterms!(model, [1, 1, 2], [1, 2, 2], [1., 1., 1.])

 # add linear constraints
add_constr!(model, [1., 2., 3.], '>', 4.)
add_constr!(model, [1., 1., 0.], '>', 1.)
update_model!(model)

optimize(model)
```


### Mixed Integer Programming

This package also supports mixed integer programming.

Problem formulation:
```
maximize x + 2 y + 5 z

s.t.  x + y + z <= 10
      x + 2 y + z <= 15
      x is continuous: 0 <= x <= 5
      y is integer: 0 <= y <= 10
      z is binary
```

##### Example 3.1: Low-level MIP API

Julia code:
```julia
using Gurobi

env = Gurobi.Env()

model = Gurobi.Model(env, "mip_01", :maximize)

 # add continuous variable
add_cvar!(model, 1., 0., 5.)  # x

 # add integer variable
add_ivar!(model, 2., 0, 10)   # y

 # add binary variable
add_bvar!(model, 5.)          # z

 # have the variables incorporated into the model
update_model!(model)

add_constr!(model, ones(3), '<', 10.)
add_constr!(model, [1., 2., 1.], '<', 15.)

optimize(model)
```

Note that you can use ``add_ivars!`` and ``add_bvars!`` to add multiple integer or binary variables in batch.

##### Example 3.2: MIP using JuMP with Gurobi

```julia
using JuMP, Gurobi

model = Model(with_optimizer(Gurobi.Optimizer))

@variables(model, begin
    0 <= x <= 5
    0 <= y <= 10, Int
    z, Bin
end)

@objective(model, Max, x + 2y + 5z)
@constraint(model, x + y + z <= 10)
@constraint(model, x + 2y + z <= 15)

optimize!(model)
```

### Quadratic constraints

The ``add_qconstr!`` function may be used to add quadratic constraints to a model.

Problem formulation:
```
maximize x + y

s.t.  x, y >= 0
      x^2 + y^2 <= 1
```

Julia code:
```julia
using Gurobi

env = Gurobi.Env()

model = Gurobi.Model(env, "qcqp_01", :maximize)

add_cvars!(model, [1., 1.], 0., Inf)
update_model!(model)

 # add_qpconstr!(model, linearindices, linearcoeffs, qrowinds, qcolinds, qcoeffs, sense, rhs)
add_qconstr!(model, [], [], [1, 2], [1, 2], [1, 1.], '<', 1.0)
update_model!(model)

optimize(model)
```

SOCP constraints of the form ``x'x <= y^2`` and ``x'x <= yz`` can be added using this method as well.

### Reusing the same Gurobi environment for multiple solves

When using this package via other packages such as [MathProgBase.jl](https://github.com/JuliaOpt/MathProgBase.jl) and [JuMP.jl](https://github.com/JuliaOpt/JuMP.jl), the default behavior is to obtain a new Gurobi license token every time a model is created and solved. If you are using Gurobi in a setting where the number of concurrent Gurobi uses is limited (e.g. ["Single-Use" or "Floating-Use" licenses](http://www.gurobi.com/products/licensing-pricing/licensing-overview)), you might instead prefer to obtain a single license token that is shared by all models that your program solves. You can do this by passing a Gurobi Environment object as the first parameter to `GurobiSolver`. For example, the follow code snippet solves multiple problems with JuMP using the same license token:

```julia
using JuMP, Gurobi
env = Gurobi.Env()

model1 = Model(with_optimizer(Gurobi.Optimizer, env))
...

# The solvers can have different options too
model2 = Model(with_optimizer(Gurobi.Optimizer, env, OutputFlag=0))
...
```
