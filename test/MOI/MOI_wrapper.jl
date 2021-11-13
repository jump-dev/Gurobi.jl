module TestMOIWrapper

using Gurobi
using Random
using Test

const MOI = Gurobi.MOI

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

const GRB_ENV = isdefined(Main, :GRB_ENV) ? Main.GRB_ENV : Gurobi.Env()

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
        exclude = String[
            # TODO(odow): investigate errors
            # Gurobi Error 10015: Cannot compute IIS on a feasible model
            # https://www.gurobi.com/documentation/9.5/refman/error_codes.html
            "test_solve_conflict_feasible",
            # SecondOrderCone does not return dual solutions. Tested below.
            "_SecondOrderCone_",
            "test_constraint_PrimalStart_DualStart_SecondOrderCone",
            "_RotatedSecondOrderCone_",
            "_GeometricMeanCone_",
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
    env = Gurobi.Env()
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

end

TestMOIWrapper.runtests()
