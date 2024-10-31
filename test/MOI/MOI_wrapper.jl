# Copyright (c) 2015 Dahua Lin, Miles Lubin, Joey Huchette, Iain Dunning, and
#   contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

module TestMOIWrapper

using Test

using Gurobi
import MathOptInterface as MOI
import Random

function runtests()
    for name in names(@__MODULE__; all = true)
        if !startswith("$(name)", "test_")
            continue
        end
        @testset "$(name)" begin
            if startswith("$(name)", "test_MULTI_ENV")
                try
                    getfield(@__MODULE__, name)()
                catch ex
                    if ex == ErrorException(
                        "Gurobi Error 10009: Failed to obtain a valid license",
                    )
                        @warn(
                            "Skipping a test because there was an issue " *
                            "creating multiple licenses. This is probably " *
                            "because you have a limited license."
                        )
                    else
                        rethrow(ex)
                    end
                end
            else
                getfield(@__MODULE__, name)()
            end
        end
    end
    return
end

const GRB_ENV =
    isdefined(Main, :GRB_ENV) ? Main.GRB_ENV : Gurobi.Env(output_flag = 0)

function test_runtests()
    model =
        MOI.Bridges.full_bridge_optimizer(Gurobi.Optimizer(GRB_ENV), Float64)
    MOI.set(model, MOI.Silent(), true)
    # We set `DualReductions = 0` so that we never return
    # `INFEASIBLE_OR_UNBOUNDED`.
    MOI.set(model, MOI.RawOptimizerAttribute("DualReductions"), 0)
    MOI.set(model, MOI.RawOptimizerAttribute("QCPDual"), 1)
    MOI.set(model, MOI.RawOptimizerAttribute("InfUnbdInfo"), 1)
    MOI.set(model, MOI.RawOptimizerAttribute("NonConvex"), 2)
    MOI.Test.runtests(
        model,
        MOI.Test.Config(atol = 1e-3, rtol = 1e-3),
        exclude = Union{String,Regex}[
            # SecondOrderCone does not return dual solutions. Tested below.
            "_SecondOrderCone_",
            "test_constraint_PrimalStart_DualStart_SecondOrderCone",
            "_RotatedSecondOrderCone_",
            "_GeometricMeanCone_",
            # Shaky tests
            "vector_nonlinear",
            "VectorNonlinearFunction",
            # Tests should be skipped due to RequirementsUnmet, but aren't
            r"^test_nonlinear_expression_hs071$",
            r"^test_nonlinear_expression_hs071_epigraph$",
            r"^test_nonlinear_expression_hs109$",
            r"^test_nonlinear_expression_hs110$",
            r"^test_nonlinear_expression_quartic$",
            r"^test_nonlinear_expression_overrides_objective$",
            r"^test_nonlinear_duals$",
        ],
    )
    MOI.Test.runtests(
        model,
        MOI.Test.Config(
            # Some of these conic ones have very low accuracy.
            atol = 2e-3,
            rtol = 1e-3,
            exclude = Any[MOI.ConstraintDual, MOI.DualObjectiveValue],
        ),
        include = String[
            "_SecondOrderCone_",
            "test_constraint_PrimalStart_DualStart_SecondOrderCone",
            "_RotatedSecondOrderCone_",
            "_GeometricMeanCone_",
        ],
    )
    return
end

function test_User_limit_handling_issue_140()
    # Verify that we return the correct status codes when a mixed-integer
    # problem has been solved to a *feasible* but not necessarily optimal
    # solution. To do that, we will set up an intentionally dumbed-down
    # Gurobi Gurobi.Optimizer (with all heuristics and pre-solve turned off) and
    # ask it to solve a classic knapsack problem. Setting SolutionLimit=1
    # forces the solver to return after its first feasible MIP solution,
    # which tests the right part of the code without relying on potentially
    # flaky or system-dependent time limits.
    m = Gurobi.Optimizer(GRB_ENV)
    MOI.set(m, MOI.Silent(), true)
    MOI.set(m, MOI.RawOptimizerAttribute("SolutionLimit"), 1)
    MOI.set(m, MOI.RawOptimizerAttribute("Heuristics"), 0)
    MOI.set(m, MOI.RawOptimizerAttribute("Presolve"), 0)
    N = 100
    x = MOI.add_variables(m, N)
    for xi in x
        MOI.add_constraint(m, xi, MOI.ZeroOne())
        MOI.set(m, MOI.VariablePrimalStart(), xi, 0.0)
    end
    # Given a collection of items with individual weights and values,
    # maximize the total value carried subject to the constraint that
    # the total weight carried is less than 10.
    Random.seed!(1)
    item_weights = rand(N)
    item_values = rand(N)
    MOI.add_constraint(
        m,
        MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(item_weights, x), 0.0),
        MOI.LessThan(10.0),
    )
    MOI.set(
        m,
        MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(),
        MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(.-item_values, x), 0.0),
    )
    MOI.optimize!(m)

    @test MOI.get(m, MOI.TerminationStatus()) == MOI.SOLUTION_LIMIT
    # We should have a primal feasible solution:
    @test MOI.get(m, MOI.PrimalStatus()) == MOI.FEASIBLE_POINT
    # But we have no dual status:
    @test MOI.get(m, MOI.DualStatus()) == MOI.NO_SOLUTION
    return
end

function test_Constant_objective_issue_111()
    m = Gurobi.Optimizer(GRB_ENV)
    x = MOI.add_variable(m)
    MOI.set(
        m,
        MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(),
        MOI.ScalarAffineFunction(MOI.ScalarAffineTerm{Float64}[], 2.0),
    )
    @test MOI.get(
        m,
        MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(),
    ).constant == 2.0
    p = Ref{Cdouble}()
    Gurobi.GRBgetdblattr(m, "ObjCon", p)
    @test p[] == 2.0

    MOI.modify(
        m,
        MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(),
        MOI.ScalarConstantChange(3.0),
    )
    @test MOI.get(
        m,
        MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(),
    ).constant == 3.0
    Gurobi.GRBgetdblattr(m.inner, "ObjCon", p)
    @test p[] == 3.0
    return
end

function test_user_provided_env()
    model_1 = Gurobi.Optimizer(GRB_ENV)
    @test model_1.env === GRB_ENV
    model_2 = Gurobi.Optimizer(GRB_ENV)
    @test model_2.env === GRB_ENV
    # Check that finalizer doesn't touch GRB_ENV when manually provided.
    finalize(model_1)
    @test GRB_ENV.ptr_env != C_NULL
    return
end

function test_MULTI_ENV()
    # Gurobi tests should pass if the function begins with test_MULTI_ENV and
    # this specific error is thrown.
    return error("Gurobi Error 10009: Failed to obtain a valid license")
end

function test_MULTI_ENV_automatic_env()
    model_1 = Gurobi.Optimizer()
    model_2 = Gurobi.Optimizer()
    @test model_1.env !== model_2.env
    # Check that env is finalized with model when not supplied manually.
    finalize(model_1)
    @test model_1.env.ptr_env == C_NULL
    return
end

function test_user_provided_env_empty()
    model = Gurobi.Optimizer(GRB_ENV)
    @test model.env === GRB_ENV
    @test GRB_ENV.ptr_env != C_NULL
    MOI.empty!(model)
    @test model.env === GRB_ENV
    @test GRB_ENV.ptr_env != C_NULL
    return
end

function test_MULTI_ENV_automatic_env_empty()
    model = Gurobi.Optimizer()
    env = model.env
    MOI.empty!(model)
    @test model.env === env
    @test env.ptr_env != C_NULL
    return
end

function test_MULTI_ENV_manual_finalize()
    env = Gurobi.Env(output_flag = 0)
    model = Gurobi.Optimizer(env)
    finalize(env)
    @test env.finalize_called
    finalize(model)
    @test env.ptr_env == C_NULL
    return
end

function test_QCPDual_1()
    model = Gurobi.Optimizer(GRB_ENV)
    MOI.set(model, MOI.Silent(), true)
    MOI.set(model, MOI.RawOptimizerAttribute("QCPDual"), 1)
    MOI.Utilities.loadfromstring!(
        model,
        """
variables: x, y, z
minobjective: 1.0 * x + 1.0 * y + 1.0 * z
c1: x + y == 2.0
c2: x + y + z >= 0.0
c3: 1.0 * x * x + -1.0 * y * y + -1.0 * z * z >= 0.0
x >= 0.0
y >= 0.0
z >= 0.0
""",
    )
    MOI.optimize!(model)
    @test MOI.get(model, MOI.TerminationStatus()) == MOI.OPTIMAL
    @test MOI.get(model, MOI.PrimalStatus()) == MOI.FEASIBLE_POINT
    @test MOI.get(model, MOI.DualStatus()) == MOI.FEASIBLE_POINT
    c1 = MOI.get(model, MOI.ConstraintIndex, "c1")
    c2 = MOI.get(model, MOI.ConstraintIndex, "c2")
    c3 = MOI.get(model, MOI.ConstraintIndex, "c3")
    @test MOI.get(model, MOI.ConstraintDual(), c1) ≈ 1.0 atol = 1e-6
    @test MOI.get(model, MOI.ConstraintDual(), c2) ≈ 0.0 atol = 1e-6
    @test MOI.get(model, MOI.ConstraintDual(), c3) ≈ 0.0 atol = 1e-6
    return
end

function test_QCPDual_default()
    model = Gurobi.Optimizer(GRB_ENV)
    MOI.set(model, MOI.Silent(), true)
    MOI.Utilities.loadfromstring!(
        model,
        """
variables: x, y, z
minobjective: 1.0 * x + 1.0 * y + 1.0 * z
c1: x + y == 2.0
c2: x + y + z >= 0.0
c3: 1.0 * x * x + -1.0 * y * y + -1.0 * z * z >= 0.0
x >= 0.0
y >= 0.0
z >= 0.0
""",
    )
    MOI.optimize!(model)
    @test MOI.get(model, MOI.TerminationStatus()) == MOI.OPTIMAL
    @test MOI.get(model, MOI.PrimalStatus()) == MOI.FEASIBLE_POINT
    @test MOI.get(model, MOI.DualStatus()) == MOI.NO_SOLUTION
    c1 = MOI.get(model, MOI.ConstraintIndex, "c1")
    c2 = MOI.get(model, MOI.ConstraintIndex, "c2")
    c3 = MOI.get(model, MOI.ConstraintIndex, "c3")
    @test_throws ErrorException MOI.get(model, MOI.ConstraintDual(), c1)
    @test_throws ErrorException MOI.get(model, MOI.ConstraintDual(), c2)
    @test_throws ErrorException MOI.get(model, MOI.ConstraintDual(), c3)
    return
end

function test_Buffered_deletion_test()
    # Check if the VarHintVal of variables (could be any attribute that is not
    # cached, the most common as: the bounds, name, and the start, are all
    # cached) are right if we remove some variables and add new ones without
    # updating the model (and that they stay correct after updating the model,
    # here done by `MOI.optimize!`, as the update needs to be done inside MOI
    # code, that is aware of the lazy updates).
    model = Gurobi.Optimizer(GRB_ENV)
    MOI.set(model, MOI.Silent(), true)
    vars = MOI.add_variables(model, 3)
    vars_obj = [7.0, 42.0, -0.5]
    obj_attr = Gurobi.VariableAttribute("Obj")
    MOI.set.(model, obj_attr, vars, vars_obj)
    @test all(MOI.is_valid.(model, vars))
    @test MOI.get.(model, obj_attr, vars) == vars_obj
    MOI.delete(model, vars[[1, 3]])
    @test !MOI.is_valid(model, vars[1])
    @test !MOI.is_valid(model, vars[3])
    fourth_var = MOI.add_variable(model)
    fourth_var_obj = -77.0
    MOI.set(model, obj_attr, fourth_var, fourth_var_obj)
    # Check before updating the model and after the delete+insert.
    @test !MOI.is_valid(model, vars[1])
    @test !MOI.is_valid(model, vars[3])
    @test MOI.is_valid(model, vars[2])
    @test vars_obj[2] == MOI.get(model, obj_attr, vars[2])
    @test MOI.is_valid(model, fourth_var)
    @test fourth_var_obj == MOI.get(model, obj_attr, fourth_var)
    # Then optimize to force an update.
    MOI.optimize!(model)
    # And check again.
    @test MOI.is_valid(model, vars[2])
    @test vars_obj[2] == MOI.get(model, obj_attr, vars[2])
    @test MOI.is_valid(model, fourth_var)
    @test fourth_var_obj == MOI.get(model, obj_attr, fourth_var)
    return
end

function test_ConstraintAttribute()
    model = Gurobi.Optimizer(GRB_ENV)
    MOI.Utilities.loadfromstring!(
        model,
        """
variables: x
minobjective: x
x >= 0.0
c: 2x >= 1.0
x in Integer()
""",
    )
    c2 = MOI.get(model, MOI.ConstraintIndex, "c")
    # Linear constraints are supported - one test for each different type.
    # Integer attribute
    MOI.set(model, Gurobi.ConstraintAttribute("Lazy"), c2, 2)
    @test MOI.get(model, Gurobi.ConstraintAttribute("Lazy"), c2) == 2
    # Real attribute
    MOI.set(model, Gurobi.ConstraintAttribute("RHS"), c2, 2.0)
    @test MOI.get(model, Gurobi.ConstraintAttribute("RHS"), c2) == 2.0
    # Char attribute
    MOI.set(model, Gurobi.ConstraintAttribute("Sense"), c2, '<')
    @test MOI.get(model, Gurobi.ConstraintAttribute("Sense"), c2) == '<'
    # String attribute
    MOI.set(model, Gurobi.ConstraintAttribute("ConstrName"), c2, "c4")
    @test MOI.get(model, Gurobi.ConstraintAttribute("ConstrName"), c2) == "c4"
    # Getting/setting a non-existing attribute.
    attr = Gurobi.ConstraintAttribute("Non-existing")
    err = ErrorException("Gurobi Error 10004: Unknown attribute 'Non-existing'")
    @test_throws err MOI.set(model, attr, c2, 1)
    @test_throws err MOI.get(model, attr, c2)
    # Setting an attribute to a value of the wrong type.
    @test_throws(
        ArgumentError(
            "Attribute Lazy requires Int64 arguments. Provided argument was of type Float64.",
        ),
        MOI.set(model, Gurobi.ConstraintAttribute("Lazy"), c2, 1.5)
    )
    return
end

function test_VariableAttribute()
    model = Gurobi.Optimizer(GRB_ENV)
    MOI.Utilities.loadfromstring!(
        model,
        """
variables: x
minobjective: x
c1: x >= 0.0
c2: 2x >= 1.0
c3: x in Integer()
""",
    )
    x = MOI.get(model, MOI.VariableIndex, "x")
    # Setting attributes of each type
    # Integer attribute
    MOI.set(model, Gurobi.VariableAttribute("VarHintPri"), x, 2)
    @test MOI.get(model, Gurobi.VariableAttribute("VarHintPri"), x) == 2
    # Real Attribute
    MOI.set(model, Gurobi.VariableAttribute("LB"), x, 2.0)
    @test MOI.get(model, Gurobi.VariableAttribute("LB"), x) == 2.0
    # Char Attribute
    MOI.set(model, Gurobi.VariableAttribute("VType"), x, 'B')
    @test MOI.get(model, Gurobi.VariableAttribute("VType"), x) == 'B'
    # String Attribute
    MOI.set(model, Gurobi.VariableAttribute("VarName"), x, "my_var")
    @test MOI.get(model, Gurobi.VariableAttribute("VarName"), x) == "my_var"
    # Things that should fail follow.
    # Getting/setting a non-existing attribute.
    attr = Gurobi.VariableAttribute("Non-existing")
    err = ErrorException("Gurobi Error 10004: Unknown attribute 'Non-existing'")
    @test_throws err MOI.set(model, attr, x, 1)
    @test_throws err MOI.get(model, attr, x)
    # Setting an attribute to a value of the wrong type.
    @test_throws(
        ArgumentError(
            "Attribute BranchPriority requires Int64 arguments. Provided argument was of type Float64.",
        ),
        MOI.set(model, Gurobi.VariableAttribute("BranchPriority"), x, 1.5)
    )
    return
end

function test_ModelAttribute()
    model = Gurobi.Optimizer(GRB_ENV)
    MOI.Utilities.loadfromstring!(
        model,
        """
variables: x
minobjective: x
c1: x >= 0.0
c2: 2x >= 1.0
c3: x in Integer()
""",
    )
    # Setting attributes of each type
    # Integer attribute
    MOI.set(model, Gurobi.ModelAttribute("ModelSense"), -1)
    @test MOI.get(model, Gurobi.ModelAttribute("ModelSense")) == -1
    # Real Attribute
    MOI.set(model, Gurobi.ModelAttribute("ObjCon"), 3.0)
    @test MOI.get(model, Gurobi.ModelAttribute("ObjCon")) == 3.0
    # String Attribute
    MOI.set(model, Gurobi.ModelAttribute("ModelName"), "My model")
    @test MOI.get(model, Gurobi.ModelAttribute("ModelName")) == "My model"
    # Things that should fail follow.
    # Getting/setting a non-existing attribute.
    attr = Gurobi.ModelAttribute("Non-existing")
    err = ErrorException("Gurobi Error 10004: Unknown attribute 'Non-existing'")
    @test_throws err MOI.set(model, attr, 1)
    @test_throws err MOI.get(model, attr)
    # Setting an attribute to a value of the wrong type.
    @test_throws(
        ArgumentError(
            "Attribute NumStart requires $(Int) arguments. Provided argument was of type Float64.",
        ),
        MOI.set(model, Gurobi.ModelAttribute("NumStart"), 4.0)
    )
    return
end

function test_Attributes()
    model = Gurobi.Optimizer(GRB_ENV)
    MOI.set(model, MOI.Silent(), true)
    MOI.optimize!(model)
    @test MOI.get(model, MOI.NodeCount()) == 0
    @test MOI.get(model, MOI.BarrierIterations()) == 0
    @test MOI.get(model, MOI.SimplexIterations()) == 0
end

function test_GRBterminate()
    model = Gurobi.Optimizer(GRB_ENV)
    MOI.set(model, MOI.Silent(), true)
    x = MOI.add_variable(model)
    MOI.set(
        model,
        Gurobi.CallbackFunction(),
        (cb_data, cb_where) -> GRBterminate(model),
    )
    MOI.optimize!(model)
    @test MOI.get(model, MOI.TerminationStatus()) == MOI.INTERRUPTED
    return
end

"""
    test_InterruptException()

This test simulates an InterruptException being thrown. It is a little
complicated due to the delayed handling of GRBterminate, which _schedules_ a
request for termination, rather than terminating immediately. This means Gurobi
may continue to call the callback after the interruption.

First, we must ensure that InterruptException() is only thrown once. Double
interrupting would interrupt our handling of the first interrupt!

Second, if the model is too simplisitic, Gurobi may be able to prove optimality
after we have interrupted, but before it has decided to actually exit the solve.
"""
function test_InterruptException()
    model = Gurobi.Optimizer(GRB_ENV)
    MOI.set(model, MOI.Silent(), true)
    x = MOI.add_variable(model)
    MOI.add_constraint(model, x, MOI.Integer())
    MOI.set(model, MOI.ObjectiveSense(), MOI.MIN_SENSE)
    MOI.set(model, MOI.ObjectiveFunction{MOI.VariableIndex}(), x)
    MOI.set(model, MOI.RawOptimizerAttribute("LazyConstraints"), 1)
    interrupt_thrown = false
    i = 0.0
    MOI.set(
        model,
        Gurobi.CallbackFunction(),
        (cb_data, cb_where) -> begin
            if cb_where != Gurobi.GRB_CB_MIPSOL
                return
            end
            MOI.submit(
                model,
                MOI.LazyConstraint(cb_data),
                MOI.ScalarAffineFunction{Float64}(
                    [MOI.ScalarAffineTerm(1.0, x)],
                    0.0,
                ),
                MOI.GreaterThan{Float64}(i),
            )
            i += 1
            if !interrupt_thrown
                interrupt_thrown = true
                throw(InterruptException())
            end
        end,
    )
    MOI.optimize!(model)
    @test MOI.get(model, MOI.TerminationStatus()) == MOI.INTERRUPTED
    return
end

function _build_basis_model()
    T = Float64
    model = Gurobi.Optimizer(GRB_ENV)
    MOI.set(model, MOI.Silent(), true)
    # Min -x
    # s.t. x + y <= 1
    # x, y >= 0

    x = MOI.add_variable(model)
    y = MOI.add_variable(model)

    cf = MOI.ScalarAffineFunction{T}(
        MOI.ScalarAffineTerm{T}.([one(T), one(T)], [x, y]),
        zero(T),
    )
    c = MOI.add_constraint(model, cf, MOI.LessThan(one(T)))

    vc1 = MOI.add_constraint(model, x, MOI.GreaterThan(zero(T)))
    vc2 = MOI.add_constraint(model, y, MOI.GreaterThan(zero(T)))
    objf = MOI.ScalarAffineFunction{T}(
        MOI.ScalarAffineTerm{T}.([-one(T), zero(T)], [x, y]),
        zero(T),
    )
    MOI.set(model, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{T}}(), objf)
    MOI.set(model, MOI.ObjectiveSense(), MOI.MIN_SENSE)
    return model, x, y, c
end

function test_set_basis()
    # The following is an indirect, brittle way to test the desired behavior.
    # Gurobi appears to not allow you to query VBasis/CBasis attributes without
    # optimizing the problem
    # (https://support.gurobi.com/hc/en-us/community/posts/360075729911-Getting-VBasis-CBasis-attributes-without-optimizing).
    # The following just verifies that the problem is solved with 0 simplex iterations
    # when we seed the optimal solution and 1 simplex iteration when we seed a
    # suboptimal feasible basis. Ideally, we would remove the call to MOI.optimize!
    # and just verify that, after setting the basis, we can get the same basis
    # statuses back.
    let
        model, x, y, c = _build_basis_model()
        MOI.set(model, MOI.VariableBasisStatus(), x, MOI.BASIC)
        MOI.set(model, MOI.VariableBasisStatus(), y, MOI.NONBASIC_AT_LOWER)
        MOI.set(model, MOI.ConstraintBasisStatus(), c, MOI.NONBASIC)
        MOI.optimize!(model)
        @test MOI.get(model, MOI.SimplexIterations()) == 0
    end
    let
        model, x, y, c = _build_basis_model()
        MOI.set(model, MOI.VariableBasisStatus(), x, MOI.NONBASIC_AT_LOWER)
        MOI.set(model, MOI.VariableBasisStatus(), y, MOI.NONBASIC_AT_LOWER)
        MOI.set(model, MOI.ConstraintBasisStatus(), c, MOI.BASIC)
        MOI.optimize!(model)
        @test MOI.get(model, MOI.SimplexIterations()) == 1
    end
    return
end

function test_add_constrained_variables()
    model = Gurobi.Optimizer(GRB_ENV)
    MOI.set(model, MOI.Silent(), true)
    set = MOI.Interval{Float64}(-1.2, 3.4)
    vi, ci = MOI.add_constrained_variable(model, set)
    @test MOI.get(model, MOI.NumberOfVariables()) == 1
    @test MOI.get(model, MOI.ListOfConstraintTypesPresent()) ==
          [(MOI.VariableIndex, MOI.Interval{Float64})]
    @test MOI.get(model, MOI.ConstraintFunction(), ci) == vi
    @test MOI.get(model, MOI.ConstraintSet(), ci) == set
    return
end

function _is_binary(x; atol = 1e-6)
    return isapprox(x, 0; atol = atol) || isapprox(x, 1; atol = atol)
end

function test_multiple_solutions()
    model = Gurobi.Optimizer(GRB_ENV)
    MOI.set(model, MOI.RawOptimizerAttribute("OutputFlag"), 0)
    MOI.set(model, MOI.RawOptimizerAttribute("Cuts"), 0)
    MOI.set(model, MOI.RawOptimizerAttribute("Presolve"), 0)
    MOI.set(model, MOI.RawOptimizerAttribute("PreCrush"), 0)
    MOI.set(model, MOI.RawOptimizerAttribute("Heuristics"), 0)
    N = 30
    x = MOI.add_variables(model, N)
    MOI.add_constraints(model, x, MOI.ZeroOne())
    MOI.set.(model, MOI.VariablePrimalStart(), x, 0.0)
    Random.seed!(1)
    item_weights, item_values = rand(N), rand(N)
    MOI.add_constraint(
        model,
        MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(item_weights, x), 0.0),
        MOI.LessThan(10.0),
    )
    MOI.set(
        model,
        MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(),
        MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(item_values, x), 0.0),
    )
    MOI.set(model, MOI.ObjectiveSense(), MOI.MAX_SENSE)
    MOI.optimize!(model)
    RC = MOI.get(model, MOI.ResultCount())
    @test RC > 1
    for n in [0, RC + 1]
        @test_throws(
            MOI.ResultIndexBoundsError,
            MOI.get(model, MOI.VariablePrimal(n), x),
        )
        @test_throws(
            MOI.ResultIndexBoundsError,
            MOI.get(model, MOI.ObjectiveValue(n)),
        )
        @test MOI.get(model, MOI.PrimalStatus(n)) == MOI.NO_SOLUTION
        @test MOI.get(model, MOI.DualStatus(n)) == MOI.NO_SOLUTION
    end
    for n in 1:RC
        xn = MOI.get(model, MOI.VariablePrimal(n), x)
        @test all(_is_binary, xn)
        @test isapprox(
            MOI.get(model, MOI.ObjectiveValue(n)),
            item_values' * xn,
            atol = 1e-6,
        )
        @test MOI.get(model, MOI.PrimalStatus(n)) == MOI.FEASIBLE_POINT
        @test MOI.get(model, MOI.DualStatus(n)) == MOI.NO_SOLUTION
    end
    return
end

function test_log_file()
    model = Gurobi.Optimizer()
    @test MOI.get(model, MOI.RawOptimizerAttribute("LogFile")) == ""
    MOI.set(model, MOI.RawOptimizerAttribute("LogFile"), "test.log")
    @test MOI.get(model, MOI.RawOptimizerAttribute("LogFile")) == "test.log"
    MOI.set(model, MOI.RawOptimizerAttribute("LogFile"), "")
    rm("test.log")
    return
end

function test_multiple_modifications()
    model = Gurobi.Optimizer(GRB_ENV)

    x = MOI.add_variables(model, 3)

    saf = MOI.ScalarAffineFunction(
        [
            MOI.ScalarAffineTerm(1.0, x[1]),
            MOI.ScalarAffineTerm(1.0, x[2]),
            MOI.ScalarAffineTerm(1.0, x[3]),
        ],
        0.0,
    )
    ci1 = MOI.add_constraint(model, saf, MOI.LessThan(1.0))
    ci2 = MOI.add_constraint(model, saf, MOI.LessThan(2.0))

    MOI.set(
        model,
        MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(),
        saf,
    )

    fc1 = MOI.get(model, MOI.ConstraintFunction(), ci1)
    @test MOI.coefficient.(fc1.terms) == [1.0, 1.0, 1.0]
    fc2 = MOI.get(model, MOI.ConstraintFunction(), ci2)
    @test MOI.coefficient.(fc2.terms) == [1.0, 1.0, 1.0]
    obj = MOI.get(
        model,
        MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(),
    )
    @test MOI.coefficient.(obj.terms) == [1.0, 1.0, 1.0]

    changes_cis = [
        MOI.ScalarCoefficientChange(MOI.VariableIndex(1), 4.0)
        MOI.ScalarCoefficientChange(MOI.VariableIndex(1), 0.5)
        MOI.ScalarCoefficientChange(MOI.VariableIndex(3), 2.0)
    ]
    MOI.modify(model, [ci1, ci2, ci2], changes_cis)

    fc1 = MOI.get(model, MOI.ConstraintFunction(), ci1)
    @test MOI.coefficient.(fc1.terms) == [4.0, 1.0, 1.0]
    fc2 = MOI.get(model, MOI.ConstraintFunction(), ci2)
    @test MOI.coefficient.(fc2.terms) == [0.5, 1.0, 2.0]

    changes_obj = [
        MOI.ScalarCoefficientChange(MOI.VariableIndex(1), 4.0)
        MOI.ScalarCoefficientChange(MOI.VariableIndex(2), 10.0)
        MOI.ScalarCoefficientChange(MOI.VariableIndex(3), 2.0)
    ]
    MOI.modify(
        model,
        MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(),
        changes_obj,
    )

    obj = MOI.get(
        model,
        MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(),
    )
    @test MOI.coefficient.(obj.terms) == [4.0, 10.0, 2.0]
end

function test_attributes_is_set_by_optimize()
    model = MOI.Utilities.CachingOptimizer(
        MOI.Utilities.UniversalFallback(MOI.Utilities.Model{Float64}()),
        Gurobi.Optimizer(GRB_ENV),
    )
    MOI.Utilities.attach_optimizer(model)
    x = MOI.add_variable(model)
    MOI.add_constraint(model, x, MOI.GreaterThan(1.5))
    MOI.add_constraint(model, x, MOI.Integer())
    c = MOI.add_constraint(model, 1.0 * x, MOI.EqualTo(2.0))
    @test MOI.get(model, Gurobi.ModelAttribute("IsMIP")) == 1
    @test MOI.get(model, Gurobi.VariableAttribute("LB"), x) == 1.5
    @test MOI.get(model, Gurobi.ConstraintAttribute("RHS"), c) == 2.0
    return
end

function test_attributes_is_copyable()
    model = MOI.Utilities.CachingOptimizer(
        MOI.Utilities.UniversalFallback(MOI.Utilities.Model{Float64}()),
        Gurobi.Optimizer(GRB_ENV),
    )
    x = MOI.add_variables(model, 2)
    MOI.set(model, Gurobi.VariableAttribute("VarHintVal"), x[1], 5.0)
    MOI.optimize!(model)
    @test MOI.get(model, MOI.TerminationStatus()) == MOI.OPTIMAL
    @test MOI.get(model, Gurobi.VariableAttribute("VarHintVal"), x[1]) == 5.0
    @test MOI.get(model, Gurobi.VariableAttribute("VarHintVal"), x[2]) == 1e101
    return
end

function test_dual_qcp_failure()
    model = MOI.Bridges.full_bridge_optimizer(Gurobi.Optimizer(), Float64)
    MOI.set(model, MOI.RawOptimizerAttribute("QCPDual"), 1)
    p = MOI.add_variables(model, 3)
    x = MOI.add_variables(model, 6)
    MOI.add_constraint.(model, p, MOI.GreaterThan(0.0))
    MOI.add_constraint.(model, x, MOI.Interval(0.0, 1.0))
    MOI.add_constraint(model, 1.0 * p[1] + p[2] + p[3], MOI.EqualTo(300.0))
    MOI.add_constraint(model, 1.0 * x[1] + x[2] + x[3], MOI.EqualTo(1.0))
    MOI.add_constraint(model, 1.0 * x[4] + x[5] + x[6], MOI.EqualTo(1.0))
    MOI.add_constraint(
        model,
        MOI.Utilities.operate(
            vcat,
            Float64,
            120.0 - 0.6p[2] + 7.2x[2] + 7.2x[5],
            2.5x[2],
            2.5x[5],
        ),
        MOI.SecondOrderCone(3),
    )
    MOI.add_constraint(
        model,
        MOI.Utilities.operate(
            vcat,
            Float64,
            120.0 - 0.6p[3] + 7.2x[3] + 7.2x[6],
            2.5x[3],
            2.5x[6],
        ),
        MOI.SecondOrderCone(3),
    )
    MOI.add_constraint(
        model,
        MOI.Utilities.operate(
            vcat,
            Float64,
            0.6p[1] - 7.2x[1] - 7.2x[4],
            2.5x[1],
            2.5x[4],
        ),
        MOI.SecondOrderCone(3),
    )
    MOI.add_constraint(
        model,
        MOI.Utilities.operate(
            vcat,
            Float64,
            0.6p[2] - 7.2x[2] - 7.2x[5],
            2.5x[2],
            2.5x[5],
        ),
        MOI.SecondOrderCone(3),
    )
    MOI.set(model, MOI.ObjectiveSense(), MOI.MIN_SENSE)
    f =
        0.15 * (p[1] - 12.0x[1] - 12.0x[4]) * (p[1] - 12.0x[1] - 12.0x[4]) +
        0.20 * (p[2] - 12.0x[2] - 12.0x[5]) * (p[2] - 12.0x[2] - 12.0x[5]) +
        0.25 * (p[3] - 12.0x[3] - 12.0x[6]) * (p[3] - 12.0x[3] - 12.0x[6]) +
        0.9x[1] * x[1] +
        1.2x[2] * x[2] +
        1.5x[3] * x[3] +
        0.9x[4] * x[4] +
        1.2x[5] * x[5] +
        1.5x[6] * x[6] +
        (1.1p[1] + 1.2p[2] + 1.3p[3]) -
        (13.2x[1] + 14.4x[2] + 15.6x[3] + 13.2x[4] + 14.4x[5] + 15.6x[6])
    MOI.set(model, MOI.ObjectiveFunction{typeof(f)}(), f)
    MOI.optimize!(model)
    MOI.get(model, MOI.DualStatus()) == MOI.NO_SOLUTION
    return
end

function test_modify_after_delete()
    model = Gurobi.Optimizer(GRB_ENV)
    x = MOI.add_variable(model)
    c = [MOI.add_constraint(model, i * x, MOI.LessThan(i)) for i in [1.0, 2.0]]
    MOI.delete(model, c[1])
    MOI.modify(model, c[2], MOI.ScalarCoefficientChange(x, 4.0))
    @test MOI.get(model, MOI.ConstraintFunction(), c[2]) ≈ 4.0 * x
    return
end

function test_modify_after_delete_plural()
    model = Gurobi.Optimizer(GRB_ENV)
    x = MOI.add_variable(model)
    c = [MOI.add_constraint(model, i * x, MOI.LessThan(i)) for i in [1.0, 2.0]]
    MOI.delete(model, c[1])
    MOI.modify(model, [c[2]], [MOI.ScalarCoefficientChange(x, 4.0)])
    @test MOI.get(model, MOI.ConstraintFunction(), c[2]) ≈ 4.0 * x
    return
end

function test_modify_after_delete_objective()
    model = Gurobi.Optimizer(GRB_ENV)
    x = MOI.add_variables(model, 2)
    MOI.set(model, MOI.ObjectiveSense(), MOI.MIN_SENSE)
    f = 1.0 * x[1] + 2.0 * x[2]
    attr = MOI.ObjectiveFunction{typeof(f)}()
    MOI.set(model, attr, f)
    MOI.delete(model, x[1])
    MOI.modify(model, attr, MOI.ScalarCoefficientChange(x[2], 3.0))
    @test MOI.get(model, attr) ≈ 3.0 * x[2]
    return
end

function test_modify_after_delete_objective_plural()
    model = Gurobi.Optimizer(GRB_ENV)
    x = MOI.add_variables(model, 2)
    MOI.set(model, MOI.ObjectiveSense(), MOI.MIN_SENSE)
    f = 1.0 * x[1] + 2.0 * x[2]
    attr = MOI.ObjectiveFunction{typeof(f)}()
    MOI.set(model, attr, f)
    MOI.delete(model, x[1])
    MOI.modify(model, attr, [MOI.ScalarCoefficientChange(x[2], 3.0)])
    @test MOI.get(model, attr) ≈ 3.0 * x[2]
    return
end

function test_attribute_TimeLimitSec()
    model = Gurobi.Optimizer(GRB_ENV)
    @test MOI.supports(model, MOI.TimeLimitSec())
    @test MOI.get(model, MOI.TimeLimitSec()) === nothing
    MOI.set(model, MOI.TimeLimitSec(), 0.0)
    @test MOI.get(model, MOI.TimeLimitSec()) == 0.0
    MOI.set(model, MOI.TimeLimitSec(), nothing)
    @test MOI.get(model, MOI.TimeLimitSec()) === nothing
    MOI.set(model, MOI.TimeLimitSec(), 1.0)
    @test MOI.get(model, MOI.TimeLimitSec()) == 1.0
    return
end

function test_last_constraint_index()
    model = Gurobi.Optimizer(GRB_ENV)
    x = MOI.add_variable(model)
    c = MOI.add_constraint(model, 1.0 * x, MOI.GreaterThan(1.0))
    @test c.value == 1
    return
end

function test_delete_indicator()
    model = Gurobi.Optimizer(GRB_ENV)
    x = MOI.add_variable(model)
    z = MOI.add_variables(model, 3)
    MOI.add_constraint.(model, z, MOI.ZeroOne())
    c = map(1:3) do i
        return MOI.add_constraint(
            model,
            MOI.Utilities.operate(vcat, Float64, z[i], 1.0 * i * x),
            MOI.Indicator{MOI.ACTIVATE_ON_ONE}(MOI.EqualTo(1.0 * i)),
        )
    end
    f = MOI.get(model, MOI.ConstraintFunction(), c[2])
    MOI.delete(model, c[1])
    MOI.delete(model, c[3])
    g = MOI.get(model, MOI.ConstraintFunction(), c[2])
    @test isapprox(f, g)
    return
end

function test_is_primal_feasible_to_tolerance()
    model = Gurobi.Optimizer(GRB_ENV)
    MOI.Utilities.loadfromstring!(
        model,
        """
        variables: x, z
        maxobjective: 2.0x + 1000.0
        c1: x + -1_333_333.12345 * z <= 0.0
        c2: z <= 0.000001
        c3: z in ZeroOne()
        """,
    )
    MOI.optimize!(model)
    @test Gurobi._is_primal_feasible_to_tolerance(model)
    return
end

function test_is_primal_feasible_to_tolerance_infeasible()
    model = Gurobi.Optimizer(GRB_ENV)
    MOI.Utilities.loadfromstring!(
        model,
        """
        variables: x, z
        maxobjective: 2.0x + 1000.0
        x + -1.0 * z <= 0.0
        z <= 0.5
        x >= 0.5
        z in ZeroOne()
        """,
    )
    MOI.optimize!(model)
    @test !Gurobi._is_primal_feasible_to_tolerance(model)
    return
end

function test_is_primal_feasible_to_tolerance_lp()
    model = Gurobi.Optimizer(GRB_ENV)
    MOI.Utilities.loadfromstring!(
        model,
        """
        variables: x
        minobjective: x
        c1: x >= 0.0
        c2: 2x >= 1.0
        """,
    )
    MOI.optimize!(model)
    @test Gurobi._is_primal_feasible_to_tolerance(model)
    return
end

function test_primal_feasible_status()
    model = Gurobi.Optimizer(GRB_ENV)
    MOI.set(model, MOI.Silent(), true)
    MOI.set(model, MOI.RawOptimizerAttribute("Heuristics"), 0)
    MOI.set(model, MOI.RawOptimizerAttribute("Presolve"), 0)
    MOI.set(model, MOI.RawOptimizerAttribute("IterationLimit"), 0)
    N = 100
    x = MOI.add_variables(model, N)
    MOI.add_constraint.(model, x, MOI.ZeroOne())
    MOI.set.(model, MOI.VariablePrimalStart(), x, 0.0)
    Random.seed!(1)
    item_weights, item_values = rand(N), rand(N)
    MOI.add_constraint(
        model,
        MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(item_weights, x), 0.0),
        MOI.LessThan(10.0),
    )
    MOI.set(model, MOI.ObjectiveSense(), MOI.MAX_SENSE)
    MOI.set(
        model,
        MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(),
        MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(item_values, x), 0.0),
    )
    MOI.optimize!(model)
    @test MOI.get(model, MOI.PrimalStatus()) == MOI.FEASIBLE_POINT
    return
end

function test_nonlinear()
    if !Gurobi._supports_nonlinear()
        return
    end
    model = Gurobi.Optimizer(GRB_ENV)
    @test MOI.supports_constraint(
        model,
        MOI.ScalarNonlinearFunction,
        MOI.LessThan{Float64},
    )
    @test MOI.supports_constraint(
        model,
        MOI.ScalarNonlinearFunction,
        MOI.GreaterThan{Float64},
    )
    @test MOI.supports_constraint(
        model,
        MOI.ScalarNonlinearFunction,
        MOI.EqualTo{Float64},
    )
    return
end

function test_nonlinear_constraint_sin()
    if !Gurobi._supports_nonlinear()
        return
    end
    model = Gurobi.Optimizer(GRB_ENV)
    MOI.set(model, MOI.Silent(), true)
    x1 = MOI.add_variable(model)
    x2 = MOI.add_variable(model)
    g = MOI.ScalarNonlinearFunction(
        :+,
        Any[MOI.ScalarNonlinearFunction(:sin, Any[2.5*x1]), 1.0*x2],
    )
    MOI.add_constraint(model, x1, MOI.GreaterThan(-1.0))
    MOI.add_constraint(model, x1, MOI.LessThan(1.0))
    MOI.add_constraint(model, x2, MOI.GreaterThan(-1.0))
    MOI.add_constraint(model, x2, MOI.LessThan(1.0))
    c = MOI.add_constraint(model, g, MOI.EqualTo(0.0))
    MOI.optimize!(model)
    x1_val = MOI.get(model, MOI.VariablePrimal(), x1)
    x2_val = MOI.get(model, MOI.VariablePrimal(), x2)
    @test ≈(sin(2.5 * x1_val) + x2_val, 0.0; atol = 1e-6)
    @test MOI.get(model, MOI.RawStatusString()) ==
          "Model was solved to optimality (subject to tolerances), and an optimal solution is available."
    return
end

function test_nonlinear_constraint_log()
    if !Gurobi._supports_nonlinear()
        return
    end
    model = Gurobi.Optimizer(GRB_ENV)
    MOI.set(model, MOI.Silent(), true)
    x = MOI.add_variable(model)
    t = MOI.add_variable(model)
    MOI.add_constraint(model, x, MOI.LessThan(2.0))
    MOI.set(model, MOI.ObjectiveSense(), MOI.MAX_SENSE)
    f = 1.0 * t
    MOI.set(model, MOI.ObjectiveFunction{typeof(f)}(), f)
    g = MOI.ScalarNonlinearFunction(
        :-,
        Any[MOI.ScalarNonlinearFunction(:log, Any[x]), t],
    )
    c = MOI.add_constraint(model, g, MOI.GreaterThan(0.0))
    MOI.optimize!(model)
    F, S = MOI.ScalarNonlinearFunction, MOI.GreaterThan{Float64}
    @test MOI.supports_constraint(model, F, S)
    @test MOI.get(model, MOI.RawStatusString()) ==
          "Model was solved to optimality (subject to tolerances), and an optimal solution is available."
    x_val = MOI.get(model, MOI.VariablePrimal(), x)
    t_val = MOI.get(model, MOI.VariablePrimal(), t)
    @test ≈(x_val, 2.0; atol = 1e-6)
    @test ≈(t_val, log(x_val); atol = 1e-6)
    @test ≈(MOI.get(model, MOI.ObjectiveValue()), t_val; atol = 1e-6)
    @test (F, S) in MOI.get(model, MOI.ListOfConstraintTypesPresent())
    @test c in MOI.get(model, MOI.ListOfConstraintIndices{F,S}())
    return
end

function test_nonlinear_constraint_unsupported()
    if !Gurobi._supports_nonlinear()
        return
    end
    model = Gurobi.Optimizer(GRB_ENV)
    x = MOI.add_variable(model)
    f = MOI.ScalarNonlinearFunction(:foo, Any[x])
    @test_throws(
        MOI.UnsupportedNonlinearOperator(:foo),
        MOI.add_constraint(model, f, MOI.GreaterThan(0.0)),
    )
    return
end

function test_nonlinear_constraint_uminus()
    if !Gurobi._supports_nonlinear()
        return
    end
    model = Gurobi.Optimizer(GRB_ENV)
    MOI.set(model, MOI.Silent(), true)
    x = MOI.add_variable(model)
    MOI.set(model, MOI.ObjectiveSense(), MOI.MAX_SENSE)
    f = 1.0 * x
    MOI.set(model, MOI.ObjectiveFunction{typeof(f)}(), f)
    g = MOI.ScalarNonlinearFunction(:-, Any[x])
    MOI.add_constraint(model, g, MOI.GreaterThan(-2.0))
    MOI.optimize!(model)
    @test MOI.get(model, MOI.RawStatusString()) ==
          "Model was solved to optimality (subject to tolerances), and an optimal solution is available."
    @test ≈(MOI.get(model, MOI.VariablePrimal(), x), 2.0; atol = 1e-3)
    @test ≈(MOI.get(model, MOI.ObjectiveValue()), 2.0; atol = 1e-3)
    return
end

function test_nonlinear_constraint_scalar_affine_function()
    if !Gurobi._supports_nonlinear()
        return
    end
    model = Gurobi.Optimizer(GRB_ENV)
    MOI.set(model, MOI.Silent(), true)
    x1 = MOI.add_variable(model)
    x2 = MOI.add_variable(model)
    x3 = MOI.add_variable(model)
    x4 = MOI.add_variable(model)
    MOI.add_constraint(model, x1, MOI.GreaterThan(0.0))
    MOI.add_constraint(model, x2, MOI.GreaterThan(0.0))
    MOI.add_constraint(model, x3, MOI.GreaterThan(0.0))
    MOI.add_constraint(model, x4, MOI.GreaterThan(0.0))
    f = 1.0 * x1 + 2.0 * x2 + 3.0 * x3
    MOI.set(model, MOI.ObjectiveSense(), MOI.MAX_SENSE)
    MOI.set(model, MOI.ObjectiveFunction{typeof(f)}(), f)
    g = MOI.ScalarNonlinearFunction(:+, Any[1.0*x1+2.0*x2+3.0*x3+4.0*x4])
    MOI.add_constraint(model, g, MOI.LessThan(6.0))
    MOI.optimize!(model)
    @test MOI.get(model, MOI.RawStatusString()) ==
          "Model was solved to optimality (subject to tolerances), and an optimal solution is available."
    @test ≈(MOI.get(model, MOI.ObjectiveValue()), 6.0; atol = 1e-3)
    return
end

function test_nonlinear_get_constraint_by_name()
    if !Gurobi._supports_nonlinear()
        return
    end
    model = Gurobi.Optimizer(GRB_ENV)
    MOI.set(model, MOI.Silent(), true)
    x = MOI.add_variable(model)
    g = MOI.ScalarNonlinearFunction(:*, Any[x, 2.0, x])
    c = MOI.add_constraint(model, g, MOI.LessThan(3.0))
    MOI.set(model, MOI.ConstraintName(), c, "c")
    d = MOI.get(model, MOI.ConstraintIndex, "c")
    @test d == c
    return
end

function test_nonlinear_constraint_delete()
    if !Gurobi._supports_nonlinear()
        return
    end
    model = Gurobi.Optimizer(GRB_ENV)
    MOI.set(model, MOI.Silent(), true)
    x = MOI.add_variable(model)
    MOI.set(model, MOI.ObjectiveSense(), MOI.MAX_SENSE)
    f = 1.0 * x
    MOI.set(model, MOI.ObjectiveFunction{typeof(f)}(), f)
    g_bad = MOI.ScalarNonlinearFunction(:exp, Any[x])
    c_bad = MOI.add_constraint(model, g_bad, MOI.GreaterThan(20.0))
    g = MOI.ScalarNonlinearFunction(:*, Any[x, 2.0, x])
    MOI.add_constraint(model, g, MOI.LessThan(3.0))
    @test MOI.is_valid(model, c_bad)
    MOI.delete(model, c_bad)
    @test !MOI.is_valid(model, c_bad)
    MOI.optimize!(model)
    @test ≈(MOI.get(model, MOI.VariablePrimal(), x), sqrt(3 / 2); atol = 1e-3)
    @test ≈(MOI.get(model, MOI.ObjectiveValue()), sqrt(3 / 2); atol = 1e-3)
    return
end

function test_nonlinear_constraint_vector_delete()
    if !Gurobi._supports_nonlinear()
        return
    end
    model = Gurobi.Optimizer(GRB_ENV)
    MOI.set(model, MOI.Silent(), true)
    x = MOI.add_variable(model)
    MOI.set(model, MOI.ObjectiveSense(), MOI.MAX_SENSE)
    f = 1.0 * x
    MOI.set(model, MOI.ObjectiveFunction{typeof(f)}(), f)
    g_bad = MOI.ScalarNonlinearFunction(:exp, Any[x])
    c_bad = MOI.add_constraint(model, g_bad, MOI.GreaterThan(20.0))
    g = MOI.ScalarNonlinearFunction(:*, Any[x, 2.0, x])
    MOI.add_constraint(model, g, MOI.LessThan(3.0))
    @test MOI.is_valid(model, c_bad)
    MOI.delete(model, [c_bad])
    @test !MOI.is_valid(model, c_bad)
    MOI.optimize!(model)
    @test ≈(MOI.get(model, MOI.VariablePrimal(), x), sqrt(3 / 2); atol = 1e-3)
    @test ≈(MOI.get(model, MOI.ObjectiveValue()), sqrt(3 / 2); atol = 1e-3)
    return
end

function test_nonlinear_pow2()
    if !Gurobi._supports_nonlinear()
        return
    end
    model = Gurobi.Optimizer(GRB_ENV)
    MOI.set(model, MOI.Silent(), true)
    # max x + y
    # s.t sqrt(^(x, 2) + ^(y, 2)) <= 1  # Use NL POW(2) operator
    # x, y >= 0
    #
    # -> x = y = 1/sqrt(2)
    x = MOI.add_variable(model)
    y = MOI.add_variable(model)
    MOI.set(model, MOI.ObjectiveSense(), MOI.MAX_SENSE)
    f = 1.0 * x + 1.0 * y
    MOI.set(model, MOI.ObjectiveFunction{typeof(f)}(), f)
    f1 = MOI.ScalarNonlinearFunction(:^, Any[x, 2])
    f2 = MOI.ScalarNonlinearFunction(:^, Any[y, 2])
    f3 = MOI.ScalarNonlinearFunction(:+, Any[f1, f2])
    g = MOI.ScalarNonlinearFunction(:sqrt, Any[f3])
    c = MOI.add_constraint(model, g, MOI.LessThan(1.0))
    MOI.optimize!(model)
    @test ≈(MOI.get(model, MOI.ObjectiveValue()), 2 / sqrt(2); atol = 1e-3)
    return
end

function test_nonlinear_scalarquadraticfunction()
    if !Gurobi._supports_nonlinear()
        return
    end
    model = Gurobi.Optimizer(GRB_ENV)
    MOI.set(model, MOI.Silent(), true)
    # max x + y
    # s.t sqrt(*(x, x) + *(y, y)) <= 1
    # x, y >= 0
    #
    # -> x = y = 1/sqrt(2)
    x = MOI.add_variable(model)
    y = MOI.add_variable(model)
    MOI.set(model, MOI.ObjectiveSense(), MOI.MAX_SENSE)
    f = 1.0 * x + 1.0 * y
    MOI.set(model, MOI.ObjectiveFunction{typeof(f)}(), f)
    f1 = MOI.ScalarNonlinearFunction(:*, Any[x, x])
    f2 = MOI.ScalarNonlinearFunction(:*, Any[y, y])
    f3 = MOI.ScalarNonlinearFunction(:+, Any[f1, f2])
    g = MOI.ScalarNonlinearFunction(:sqrt, Any[f3])
    c = MOI.add_constraint(model, g, MOI.LessThan(1.0))
    MOI.optimize!(model)
    @test ≈(MOI.get(model, MOI.ObjectiveValue()), 2 / sqrt(2); atol = 1e-3)
    return
end

function test_ModelName_too_long()
    model = Gurobi.Optimizer(GRB_ENV)
    @test_throws(
        MOI.SetAttributeNotAllowed{MOI.Name},
        MOI.set(model, MOI.Name(), "a"^(GRB_MAX_NAMELEN + 1)),
    )
    name = "a"^GRB_MAX_NAMELEN
    MOI.set(model, MOI.Name(), name)
    @test MOI.get(model, MOI.Name()) == name
    return
end

function test_VarName_too_long()
    model = Gurobi.Optimizer(GRB_ENV)
    x = MOI.add_variable(model)
    name = "a"^256
    MOI.set(model, MOI.VariableName(), x, name)
    @test MOI.get(model, MOI.VariableName(), x) == name
    MOI.set(model, MOI.VariableName(), x, "x")
    @test MOI.get(model, MOI.VariableName(), x) == "x"
    MOI.set(model, MOI.VariableName(), x, "")
    @test MOI.get(model, MOI.VariableName(), x) == ""
    return
end

function test_ConstrName_too_long()
    model = Gurobi.Optimizer(GRB_ENV)
    x = MOI.add_variable(model)
    c = MOI.add_constraint(model, 1.0 * x, MOI.EqualTo(2.0))
    name = "a"^256
    MOI.set(model, MOI.ConstraintName(), c, name)
    @test MOI.get(model, MOI.ConstraintName(), c) == name
    MOI.set(model, MOI.ConstraintName(), c, "c")
    @test MOI.get(model, MOI.ConstraintName(), c) == "c"
    MOI.set(model, MOI.ConstraintName(), c, "")
    @test MOI.get(model, MOI.ConstraintName(), c) == ""
    return
end

end  # TestMOIWrapper

TestMOIWrapper.runtests()
