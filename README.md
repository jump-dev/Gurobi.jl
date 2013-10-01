## Gurobi.jl

The [Gurobi](http://www.gurobi.com) Optimizer is a commercial optimization solver for a variety of mathematical programming problems, including linear programming (LP), quadratic programming (QP), quadratically constrained programming (QCP), mixed integer linear programming (MILP), mixed-integer quadratic programming (MIQP), and mixed-integer quadratically constrained programming (MIQCP).

The Gurobi solver is considered one of the best solvers (in terms of performance and success rate of tackling hard problems) in math programming, and its performance is comparable to (and sometimes superior to) CPLEX.
While in general it would be expensive to purchase a Gurobi license, academic users can get a license for free. 

This package is a wrapper of the Gurobi solver (through its C interface). Currently, this package supports the following types of problems:

* Linear programming (LP)
* Mixed Integer Linear Programming (MILP)
* Quadratic programming (QP)
* Mixed Integer Quadratic Programming (MIQP)
* Quadratically constrained quadratic programming (QCQP)
* Second order cone programming (SOCP)
* Mixed integer second order cone programming (MISOCP)

### Installation

Here is the procedure to setup this package:

1. Obtain a license of Gurobi and install Gurobi solver, following the instructions on [Gurobi's website](http://www.gurobi.com).

2. Install this package using ``Pkg.add("Gurobi")``.

3. Make sure the ``GUROBI_HOME`` environmental variable is set to the path of the Gurobi directory. This is part of a standard installation. The Gurobi library will be searched for in ``GUROBI_HOME/lib``. If the library is not found, check that your version is listed in ``deps/build.jl``.

4. Now, you can start using it.


### Examples

The usage of this package is straight forward. Here, it demonstrates the use of this package through several examples.

#### Example 1: Low-level Linear Programming

Problem formulation:
```
maximize x + y

s.t. 50 x + 24 y <= 2400
     30 x + 33 y <= 2100
     x >= 45, y >= 5
```

Julia code:
```julia
using Gurobi

 # creates an environment, which captures the solver setting 
env = Gurobi.Env()

 # creates an empty model ("lp_01" is the model name)
model = gurobi_model(env, "lp_01", :maximize)

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

 # show results
sol = get_solution(model)
println("soln = $(sol)")

objv = get_objval(model)
println("objv = $(objv)")
```

You may also add variables and constraints in batch, as:

```julia

 # add mutliple variables in batch
 # add_cvar!(model, obj_coefs, lower_bound, upper_bound)
add_cvars!(model, [1., 1.], [45., 5.], nothing)

 # add multiple constraints in batch 
 # add_constrs!(model, rowbegin, cind, coeffs, sense, rhs)
 # add_constrs!(model, Cint[1, 3], Cint[1, 2, 1, 2], 
 #    [50., 24., 30., 33.], '<', [2400., 2100.])
```

In ``add_constrs``, the left hand side matrix should be input in the form of compressed sparse rows. Specifically, the begining index of the i-th row is given by ``rowbegin[i]``. In the example above, we have

```julia
rowbegin = Cint[1, 3]
cind = Cint[1, 2, 1, 2] 
coeffs = [50., 24., 30., 33.]
```

Here, ``Cint[1, 3]`` implies the range ``1:2`` is for the first row. Therefore, the first row have nonzeros at positions ``[1, 2]`` and their values are ``[50., 24.]``. Likewise, the second row have nonzeros at ``[1, 2]`` and their values are ``[30., 33.]``.

#### Example 2: Linear programming (MATLAB-like style)

You may also specify the entire problem in one function call, using the 
solver-independent **[MathProgBase]** package. See that package
for more information.

[MathProgBase]: https://github.com/mlubin/MathProgBase.jl

Problem formulation:
```
maximize 1000 x + 350 y

s.t. x >= 30, y >= 0
     -1. x + 1.5 y <= 0
     12. x + 8.  y <= 1000
     1000 x + 300 y <= 70000
```

Julia code:
```julia
using MathProgBase

f = [1000., 350.]
A = [-1. 1.5; 12. 8.; 1000. 300.]
b = [0., 1000., 70000.]
lb = [0., 30.]
 
# linprog always minimizes, so we flip the objective
solution = linprog(-f,A,'<',b,lb,Inf, LPSolver(:Gurobi))
```

#### Example 3: Linear programming (Algebraic model)

Using **[JuMP]**, we can specify linear programming problems using a more
natural algebraic approach.

[JuMP]: https://github.com/IainNZ/JuMP.jl

Problem formulation:
```
maximize 1000 x + 350 y

s.t. x >= 30, y >= 0
     -1. x + 1.5 y <= 0
     12. x + 8.  y <= 1000
     1000 x + 300 y <= 70000
```

Julia code:
```julia
using JuMP

m = Model(:Max,lpsolver=LPSolver(:Gurobi))

@defVar(m, x >= 30)
@defVar(m, y >= 0)

@setObjective(m, 1000x + 350y)
@addConstraint(m, -x + 1.5y <= 0)
@addConstraint(m, 12x + 8y <= 1000)
@addConstraint(m, 1000x + 300y <= 70000)

status = solve(m)
println("Optimal objective: ",getObjectiveValue(m), 
	". x = ", getValue(x), " y = ", getValue(y))
```


#### Example 4: Low-level Quadratic programming  

To construct a QP model using the low-level interface, you have to add QP terms to a model using ``add_qpterms``

Problem formulation:
```
minimize 2 x^2 + y^2 + xy + x + y

s.t.  x, y >= 0
      x + y = 1
```

Julia code:
```julia
env = Gurobi.Env()

model = gurobi_model(env, "qp_02")

add_cvars!(model, [1., 1.], 0., Inf)
update_model!(model)

 # add quadratic terms: 2 x^2, x * y, y^2
 # add_qpterms!(model, rowinds, colinds, coeffs)
add_qpterms!(model, [1, 1, 2], [1, 2, 2], [2., 1., 1.])
add_constr!(model, [1., 1.], '=', 1.)
update_model!(model)

optimize(model)
```

#### Example 5: Quadratic programming (MATLAB-like style)

As MathProgBase does not yet support quadratic programming,
this package provides a ``qp_model`` function to construct QP problems in a style like MATLAB's ``quadprog``.

Problem formulation:
```
minimize x^2 + xy + y^2 + yz + z^2

s.t.  x + 2 y + 3 z >= 4
      x +   y       >= 1
```

Julia code:
```julia
env = Gurobi.Env()

H = [2. 1. 0.; 1. 2. 1.; 0. 1. 2.]
f = zeros(3)
A = -[1. 2. 3.; 1. 1. 0.]
b = -[4., 1.]

model = qp_model(env, "qp_02", H, f, A, b)
optimize(model)
```

#### Example 6: Mixed Integer Programming

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

Julia code:
```julia
env = Gurobi.Env()
model = gurobi_model(env, "mip_01", :maximize)

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

The same model can be more easily expressed by using JuMP:

```julia
using JuMP

m = Model(:Max, mipsolver=MIPSolver(:Gurobi))

@defVar(m, 0 <= x <= 5)
@defVar(m, 0 <= y <= 10, Int)
@defVar(m, z, Bin)

@setObjective(m, x + 2y + 5z)
@addConstraint(m, x + y + z <= 10)
@addConstraint(m, x + 2y + z <= 15)

solve(m)
```

#### Example 7: Quadratic constraints

The ``add_qconstr!`` function may be used to add quadratic constraints to a model.

Problem formulation:
```
maximize x + y

s.t.  x, y >= 0
      x^2 + y^2 <= 1
```

Julia code:
```julia
env = Gurobi.Env()

model = gurobi_model(env, "qcqp_01", :maximize)

add_cvars!(model, [1., 1.], 0., Inf)
update_model!(model)

 # add_qpconstr!(model, linearindices, linearcoeffs, qrowinds, qcolinds, qcoeffs, sense, rhs)
add_qconstr!(model, [], [], [1, 2], [1, 2], [1, 1.], '<', 1.0)
update_model!(model)

optimize(model)
```

SOCP constraints of the form ``x'x <= y^2`` and ``x'x <= yz`` can be added using this method as well.

### Parameter Settings

In Gurobi, solver parameters are encapsulated by the ``Env`` instance. This package provides functions to get and set parameters

```julia
get_int_param(env, name)      # get an integer parameter
get_dbl_param(env, name)      # get a real-valued parameter

set_int_param!(env, name, v)   # set an integer parameter
set_dbl_param!(env, name, v)   # set a real-valued parameter
```

You may refer to Gurobi's [Parameter Reference](http://www.gurobi.com/documentation/5.0/reference-manual/node653) for the whole list of parameters. 

Here are some simple examples
```julia
set_int_param!(env, "Method", 2)   # choose to use Barrier method
set_dbl_param!(env, "IterationLimit", 100.) # set the maximum iterations (for Simplex)
```

These parameters may be used directly with the Gurobi ``LPSolver`` and ``MIPSolver``
objects from MathProgBase. For example:
```julia
solver = LPSolver(:Gurobi, Method=2)
solver = LPSolver(:Gurobi, IterationLimit=100.)
```

Note that type of the value of the parameter is used to infer whether it corresponds to
an integer or real-valued parameter in Gurobi. ``LPSolver(:Gurobi, IterationLimit=100)``
will therefore cause an error, because ``100`` is an integer and ``IterationLimit``
is a real-valued parameter.

