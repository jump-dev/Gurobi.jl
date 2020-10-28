module TestMOIWrapper

using Gurobi
using Random
using Test

const MOI  = Gurobi.MOI
const MOIT = MOI.Test

const GRB_ENV = isdefined(Main, :GRB_ENV) ? Main.GRB_ENV : Gurobi.Env()

const OPTIMIZER = MOI.Bridges.full_bridge_optimizer(
    begin
        model = Gurobi.Optimizer(GRB_ENV)
        MOI.set(model, MOI.Silent(), true)
        # We set `DualReductions = 0` so that we never return
        # `INFEASIBLE_OR_UNBOUNDED`.
        MOI.set(model, MOI.RawParameter("DualReductions"), 0)
        MOI.set(model, MOI.RawParameter("QCPDual"), 1)
        MOI.set(model, MOI.RawParameter("InfUnbdInfo"), 1)
        MOI.set(model, MOI.RawParameter("NonConvex"), 2)
        model
    end,
    Float64
)
const CONFIG = MOIT.TestConfig()

function test_basic_constraint_tests()
    MOIT.basic_constraint_tests(OPTIMIZER, CONFIG; exclude = [
        (MOI.VectorOfVariables, MOI.SecondOrderCone),
        (MOI.VectorOfVariables, MOI.RotatedSecondOrderCone),
        (MOI.VectorOfVariables, MOI.GeometricMeanCone),
        (MOI.VectorAffineFunction{Float64}, MOI.SecondOrderCone),
        (MOI.VectorAffineFunction{Float64}, MOI.RotatedSecondOrderCone),
        (MOI.VectorAffineFunction{Float64}, MOI.GeometricMeanCone),
        (MOI.VectorQuadraticFunction{Float64}, MOI.SecondOrderCone),
        (MOI.VectorQuadraticFunction{Float64}, MOI.RotatedSecondOrderCone),
        (MOI.VectorQuadraticFunction{Float64}, MOI.GeometricMeanCone),
    ])
    # TODO(odow): bugs deleting SOC variables. See also the
    # `delete_soc_variables` test.
    MOIT.basic_constraint_tests(
        OPTIMIZER,
        CONFIG;
        include = [
            (MOI.VectorOfVariables, MOI.SecondOrderCone),
            (MOI.VectorOfVariables, MOI.RotatedSecondOrderCone),
            (MOI.VectorOfVariables, MOI.GeometricMeanCone),
            (MOI.VectorAffineFunction{Float64}, MOI.SecondOrderCone),
            (MOI.VectorAffineFunction{Float64}, MOI.RotatedSecondOrderCone),
            (MOI.VectorAffineFunction{Float64}, MOI.GeometricMeanCone),
            (MOI.VectorQuadraticFunction{Float64}, MOI.SecondOrderCone),
            (MOI.VectorQuadraticFunction{Float64}, MOI.RotatedSecondOrderCone),
            (MOI.VectorQuadraticFunction{Float64}, MOI.GeometricMeanCone),
        ],
        delete = false
    )
end

function test_unittest()
    MOIT.unittest(OPTIMIZER, MOIT.TestConfig(atol=1e-6), [
        # TODO(odow): bug! We can't delete a vector of variables if one is in
        # a second order cone.
        "delete_soc_variables",
    ])
end

function test_modificationtest()
    MOIT.modificationtest(OPTIMIZER, CONFIG)
end

function test_contlineartest()
    MOIT.contlineartest(OPTIMIZER, MOIT.TestConfig(basis = true))
end

function test_contquadratictest()
    MOIT.contquadratictest(OPTIMIZER, MOIT.TestConfig(atol=1e-3, rtol=1e-3))
end

function test_conictest()
    MOIT.lintest(OPTIMIZER, CONFIG)
    MOIT.soctest(OPTIMIZER, MOIT.TestConfig(duals = false, atol=1e-3), ["soc3"])
    MOIT.soc3test(
        OPTIMIZER,
        MOIT.TestConfig(duals = false, infeas_certificates = false, atol = 1e-3)
    )
    MOIT.rsoctest(OPTIMIZER, MOIT.TestConfig(duals = false, atol=5e-3))
    MOIT.geomeantest(OPTIMIZER, MOIT.TestConfig(duals = false, atol=1e-3))
end

function test_intlinear()
    MOIT.intlineartest(OPTIMIZER, CONFIG)
end

function test_solvername()
    @test MOI.get(OPTIMIZER, MOI.SolverName()) == "Gurobi"
end

function test_default_objective_test()
    MOIT.default_objective_test(OPTIMIZER)
end

function test_default_status_test()
    MOIT.default_status_test(OPTIMIZER)
end

function test_nametest()
    MOIT.nametest(OPTIMIZER)
end

function test_validtest()
    MOIT.validtest(OPTIMIZER)
end

function test_emptytest()
    MOIT.emptytest(OPTIMIZER)
end

function test_orderedindicestest()
    MOIT.orderedindicestest(OPTIMIZER)
end

function test_copytest()
    MOIT.copytest(
        OPTIMIZER,
        MOI.Bridges.full_bridge_optimizer(Gurobi.Optimizer(GRB_ENV), Float64)
    )
end

function test_scalar_function_constant_not_zero()
    MOIT.scalar_function_constant_not_zero(OPTIMIZER)
end

function test_start_values_test()
    model = Gurobi.Optimizer(GRB_ENV)
    MOI.set(model, MOI.Silent(), true)
    x = MOI.add_variables(model, 2)
    @test MOI.supports(model, MOI.VariablePrimalStart(), MOI.VariableIndex)
    @test MOI.get(model, MOI.VariablePrimalStart(), x[1]) === nothing
    @test MOI.get(model, MOI.VariablePrimalStart(), x[2]) === nothing
    Gurobi._update_if_necessary(model)
    p = Ref{Cdouble}()
    Gurobi.GRBgetdblattrelement(model.inner, "Start", Cint(0), p)
    @test p[] == Gurobi.GRB_UNDEFINED
    MOI.set(model, MOI.VariablePrimalStart(), x[1], 1.0)
    MOI.set(model, MOI.VariablePrimalStart(), x[2], nothing)
    @test MOI.get(model, MOI.VariablePrimalStart(), x[1]) == 1.0
    @test MOI.get(model, MOI.VariablePrimalStart(), x[2]) === nothing
    Gurobi._update_if_necessary(model)
    Gurobi.GRBgetdblattrelement(model.inner, "Start", Cint(1), p)
    @test p[] == Gurobi.GRB_UNDEFINED
    MOI.optimize!(model)
    @test MOI.get(model, MOI.ObjectiveValue()) == 0.0
    # We don't support ConstraintDualStart or ConstraintPrimalStart yet.
    # @test_broken MOIT.start_values_test(Gurobi.Optimizer(GRB_ENV), OPTIMIZER)
end

function test_supports_constrainttest()
    # supports_constrainttest needs VectorOfVariables-in-Zeros,
    # MOIT.supports_constrainttest(Gurobi.Optimizer(GRB_ENV), Float64, Float32)
    # but supports_constrainttest is broken via bridges:
    MOI.empty!(OPTIMIZER)
    MOI.add_variable(OPTIMIZER)
    @test  MOI.supports_constraint(OPTIMIZER, MOI.SingleVariable, MOI.EqualTo{Float64})
    @test  MOI.supports_constraint(OPTIMIZER, MOI.ScalarAffineFunction{Float64}, MOI.EqualTo{Float64})
    # This test is broken for some reason:
    @test_broken !MOI.supports_constraint(OPTIMIZER, MOI.ScalarAffineFunction{Int}, MOI.EqualTo{Float64})
    @test !MOI.supports_constraint(OPTIMIZER, MOI.ScalarAffineFunction{Int}, MOI.EqualTo{Int})
    @test !MOI.supports_constraint(OPTIMIZER, MOI.SingleVariable, MOI.EqualTo{Int})
    @test  MOI.supports_constraint(OPTIMIZER, MOI.VectorOfVariables, MOI.Zeros)
    @test !MOI.supports_constraint(OPTIMIZER, MOI.VectorOfVariables, MOI.EqualTo{Float64})
    @test !MOI.supports_constraint(OPTIMIZER, MOI.SingleVariable, MOI.Zeros)
    @test !MOI.supports_constraint(OPTIMIZER, MOI.VectorOfVariables, MOIT.UnknownVectorSet)
end

function test_set_lower_bound_twice()
    MOIT.set_lower_bound_twice(OPTIMIZER, Float64)
end

function test_set_upper_bound_twice()
    MOIT.set_upper_bound_twice(OPTIMIZER, Float64)
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
    MOI.set(m, MOI.RawParameter("SolutionLimit"), 1)
    MOI.set(m, MOI.RawParameter("Heuristics"), 0)
    MOI.set(m, MOI.RawParameter("Presolve"), 0)
    N = 100
    x = MOI.add_variables(m, N)
    for xi in x
        MOI.add_constraint(m, MOI.SingleVariable(xi), MOI.ZeroOne())
        MOI.set(m, MOI.VariablePrimalStart(), xi, 0.0)
    end
    # Given a collection of items with individual weights and values,
    # maximize the total value carried subject to the constraint that
    # the total weight carried is less than 10.
    Random.seed!(1)
    item_weights = rand(N)
    item_values = rand(N)
    MOI.add_constraint(m,
        MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(item_weights, x), 0.0),
        MOI.LessThan(10.0))
    MOI.set(m, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(),
        MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(.-item_values, x), 0.0))
    MOI.optimize!(m)

    @test MOI.get(m, MOI.TerminationStatus()) == MOI.SOLUTION_LIMIT
    # We should have a primal feasible solution:
    @test MOI.get(m, MOI.PrimalStatus()) == MOI.FEASIBLE_POINT
    # But we have no dual status:
    @test MOI.get(m, MOI.DualStatus()) == MOI.NO_SOLUTION
end

function test_Constant_objective_issue_111()
    m = Gurobi.Optimizer(GRB_ENV)
    x = MOI.add_variable(m)
    MOI.set(m, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(),
        MOI.ScalarAffineFunction(MOI.ScalarAffineTerm{Float64}[], 2.0))
    @test MOI.get(m, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}()).constant == 2.0
    p = Ref{Cdouble}()
    Gurobi.GRBgetdblattr(m, "ObjCon", p)
    @test p[] == 2.0

    MOI.modify(m, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(), MOI.ScalarConstantChange(3.0))
    @test MOI.get(m, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}()).constant == 3.0
    Gurobi.GRBgetdblattr(m.inner, "ObjCon", p)
    @test p[] == 3.0
end

function test_user_provided_env()
    model_1 = Gurobi.Optimizer(GRB_ENV)
    @test model_1.env === GRB_ENV
    model_2 = Gurobi.Optimizer(GRB_ENV)
    @test model_2.env === GRB_ENV
    # Check that finalizer doesn't touch GRB_ENV when manually provided.
    finalize(model_1)
    @test GRB_ENV.ptr_env != C_NULL
end

function test_MULTI_ENV()
    # Gurobi tests should pass if the function begins with test_MULTI_ENV and
    # this specific error is thrown.
    error("Gurobi Error 10009: Failed to obtain a valid license")
end

function test_MULTI_ENV_automatic_env()
    model_1 = Gurobi.Optimizer()
    model_2 = Gurobi.Optimizer()
    @test model_1.env !== model_2.env
    # Check that env is finalized with model when not supplied manually.
    finalize(model_1)
    @test model_1.env.ptr_env == C_NULL
end

function test_user_provided_env_empty()
    model = Gurobi.Optimizer(GRB_ENV)
    @test model.env === GRB_ENV
    @test GRB_ENV.ptr_env != C_NULL
    MOI.empty!(model)
    @test model.env === GRB_ENV
    @test GRB_ENV.ptr_env != C_NULL
end

function test_MULTI_ENV_automatic_env_empty()
    model = Gurobi.Optimizer()
    env = model.env
    MOI.empty!(model)
    @test model.env === env
    @test env.ptr_env != C_NULL
end

function test_MULTI_ENV_manual_finalize()
    env = Gurobi.Env()
    model = Gurobi.Optimizer(env)
    finalize(env)
    @test env.finalize_called
    finalize(model)
    @test env.ptr_env == C_NULL
end

function test_Conflict_refiner_bound_bound()
    model = Gurobi.Optimizer(GRB_ENV)
    MOI.set(model, MOI.Silent(), true)
    x = MOI.add_variable(model)
    c1 = MOI.add_constraint(model, MOI.SingleVariable(x), MOI.GreaterThan(2.0))
    c2 = MOI.add_constraint(model, MOI.SingleVariable(x), MOI.LessThan(1.0))

    # Getting the results before the conflict refiner has been called must return an error.
    @test MOI.get(model, Gurobi.ConflictStatus()) == -1
    @test MOI.get(model, MOI.ConflictStatus()) == MOI.COMPUTE_CONFLICT_NOT_CALLED
    @test_throws ErrorException MOI.get(model, MOI.ConstraintConflictStatus(), c1)

    # Once it's called, no problem.
    MOI.compute_conflict!(model)
    @test MOI.get(model, Gurobi.ConflictStatus()) == 0
    @test MOI.get(model, MOI.ConflictStatus()) == MOI.CONFLICT_FOUND
    @test MOI.get(model, MOI.ConstraintConflictStatus(), c1) == MOI.IN_CONFLICT
    @test MOI.get(model, MOI.ConstraintConflictStatus(), c2) == MOI.IN_CONFLICT
end

function test_Conflict_refiner_bound_affine()
    model = Gurobi.Optimizer(GRB_ENV)
    MOI.set(model, MOI.Silent(), true)
    x = MOI.add_variable(model)
    c1 = MOI.add_constraint(model, MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.([1.0], [x]), 0.0), MOI.GreaterThan(2.0))
    c2 = MOI.add_constraint(model, MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.([1.0], [x]), 0.0), MOI.LessThan(1.0))

    # Getting the results before the conflict refiner has been called must return an error.
    @test MOI.get(model, Gurobi.ConflictStatus()) == -1
    @test MOI.get(model, MOI.ConflictStatus()) == MOI.COMPUTE_CONFLICT_NOT_CALLED
    @test_throws ErrorException MOI.get(model, MOI.ConstraintConflictStatus(), c1)

    # Once it's called, no problem.
    MOI.compute_conflict!(model)
    @test MOI.get(model, Gurobi.ConflictStatus()) == 0
    @test MOI.get(model, MOI.ConflictStatus()) == MOI.CONFLICT_FOUND
    @test MOI.get(model, MOI.ConstraintConflictStatus(), c1) == MOI.IN_CONFLICT
    @test MOI.get(model, MOI.ConstraintConflictStatus(), c2) == MOI.IN_CONFLICT
end

function test_Conflict_refiner_invalid_interval()
    model = Gurobi.Optimizer(GRB_ENV)
    MOI.set(model, MOI.Silent(), true)
    x = MOI.add_variable(model)
    c1 = MOI.add_constraint(
        model, MOI.SingleVariable(x), MOI.Interval(1.0, 0.0)
    )
    # Getting the results before the conflict refiner has been called must return an error.
    @test MOI.get(model, Gurobi.ConflictStatus()) == -1
    @test MOI.get(model, MOI.ConflictStatus()) == MOI.COMPUTE_CONFLICT_NOT_CALLED
    @test_throws ErrorException MOI.get(model, MOI.ConstraintConflictStatus(), c1)

    # Once it's called, no problem.
    MOI.compute_conflict!(model)
    @test MOI.get(model, Gurobi.ConflictStatus()) == 0
    @test MOI.get(model, MOI.ConflictStatus()) == MOI.CONFLICT_FOUND
    @test MOI.get(model, MOI.ConstraintConflictStatus(), c1) == MOI.IN_CONFLICT
end

function test_Conflict_refiner_affine_affine()
    model = Gurobi.Optimizer(GRB_ENV)
    MOI.set(model, MOI.Silent(), true)
    x = MOI.add_variable(model)
    y = MOI.add_variable(model)
    b1 = MOI.add_constraint(model, MOI.SingleVariable(x), MOI.GreaterThan(0.0))
    b2 = MOI.add_constraint(model, MOI.SingleVariable(y), MOI.GreaterThan(0.0))
    cf1 = MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.([1.0, 1.0], [x, y]), 0.0)
    c1 = MOI.add_constraint(model, cf1, MOI.LessThan(-1.0))
    cf2 = MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.([1.0, -1.0], [x, y]), 0.0)
    c2 = MOI.add_constraint(model, cf2, MOI.GreaterThan(1.0))

    # Getting the results before the conflict refiner has been called must return an error.
    @test MOI.get(model, Gurobi.ConflictStatus()) == -1
    @test MOI.get(model, MOI.ConflictStatus()) == MOI.COMPUTE_CONFLICT_NOT_CALLED
    @test_throws ErrorException MOI.get(model, MOI.ConstraintConflictStatus(), c1)

    # Once it's called, no problem.
    MOI.compute_conflict!(model)
    @test MOI.get(model, Gurobi.ConflictStatus()) == 0
    @test MOI.get(model, MOI.ConflictStatus()) == MOI.CONFLICT_FOUND
    @test MOI.get(model, MOI.ConstraintConflictStatus(), b1) == MOI.IN_CONFLICT
    @test MOI.get(model, MOI.ConstraintConflictStatus(), b2) == MOI.IN_CONFLICT
    @test MOI.get(model, MOI.ConstraintConflictStatus(), c1) == MOI.IN_CONFLICT
    @test MOI.get(model, MOI.ConstraintConflictStatus(), c2) == MOI.NOT_IN_CONFLICT
end

function test_Conflict_refiner_equalto()
    model = Gurobi.Optimizer(GRB_ENV)
    MOI.set(model, MOI.Silent(), true)
    x = MOI.add_variable(model)
    y = MOI.add_variable(model)
    b1 = MOI.add_constraint(model, MOI.SingleVariable(x), MOI.GreaterThan(0.0))
    b2 = MOI.add_constraint(model, MOI.SingleVariable(y), MOI.GreaterThan(0.0))
    cf1 = MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.([1.0, 1.0], [x, y]), 0.0)
    c1 = MOI.add_constraint(model, cf1, MOI.EqualTo(-1.0))
    cf2 = MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.([1.0, -1.0], [x, y]), 0.0)
    c2 = MOI.add_constraint(model, cf2, MOI.GreaterThan(1.0))

    # Getting the results before the conflict refiner has been called must return an error.
    @test MOI.get(model, Gurobi.ConflictStatus()) == -1
    @test MOI.get(model, MOI.ConflictStatus()) == MOI.COMPUTE_CONFLICT_NOT_CALLED
    @test_throws ErrorException MOI.get(model, MOI.ConstraintConflictStatus(), c1)

    # Once it's called, no problem.
    MOI.compute_conflict!(model)
    @test MOI.get(model, Gurobi.ConflictStatus()) == 0
    @test MOI.get(model, MOI.ConflictStatus()) == MOI.CONFLICT_FOUND
    @test MOI.get(model, MOI.ConstraintConflictStatus(), b1) == MOI.IN_CONFLICT
    @test MOI.get(model, MOI.ConstraintConflictStatus(), b2) == MOI.IN_CONFLICT
    @test MOI.get(model, MOI.ConstraintConflictStatus(), c1) == MOI.IN_CONFLICT
    @test MOI.get(model, MOI.ConstraintConflictStatus(), c2) == MOI.NOT_IN_CONFLICT
end

function test_Conflict_refiner_outside_conflict()
    model = Gurobi.Optimizer(GRB_ENV)
    MOI.set(model, MOI.Silent(), true)
    x = MOI.add_variable(model)
    y = MOI.add_variable(model)
    z = MOI.add_variable(model)
    b1 = MOI.add_constraint(model, MOI.SingleVariable(x), MOI.GreaterThan(0.0))
    b2 = MOI.add_constraint(model, MOI.SingleVariable(y), MOI.GreaterThan(0.0))
    b3 = MOI.add_constraint(model, MOI.SingleVariable(z), MOI.GreaterThan(0.0))
    cf1 = MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.([1.0, 1.0], [x, y]), 0.0)
    c1 = MOI.add_constraint(model, cf1, MOI.LessThan(-1.0))
    cf2 = MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.([1.0, -1.0, 1.0], [x, y, z]), 0.0)
    c2 = MOI.add_constraint(model, cf2, MOI.GreaterThan(1.0))

    # Getting the results before the conflict refiner has been called must return an error.
    @test MOI.get(model, Gurobi.ConflictStatus()) == -1
    @test MOI.get(model, MOI.ConflictStatus()) == MOI.COMPUTE_CONFLICT_NOT_CALLED
    @test_throws ErrorException MOI.get(model, MOI.ConstraintConflictStatus(), c1)

    # Once it's called, no problem.
    MOI.compute_conflict!(model)
    @test MOI.get(model, Gurobi.ConflictStatus()) == 0
    @test MOI.get(model, MOI.ConflictStatus()) == MOI.CONFLICT_FOUND
    @test MOI.get(model, MOI.ConstraintConflictStatus(), b1) == MOI.IN_CONFLICT
    @test MOI.get(model, MOI.ConstraintConflictStatus(), b2) == MOI.IN_CONFLICT
    @test MOI.get(model, MOI.ConstraintConflictStatus(), b3) == MOI.NOT_IN_CONFLICT
    @test MOI.get(model, MOI.ConstraintConflictStatus(), c1) == MOI.IN_CONFLICT
    @test MOI.get(model, MOI.ConstraintConflictStatus(), c2) == MOI.NOT_IN_CONFLICT
end

function test_Conflict_refiner_no_conflict()
    model = Gurobi.Optimizer(GRB_ENV)
    MOI.set(model, MOI.Silent(), true)
    x = MOI.add_variable(model)
    c1 = MOI.add_constraint(model, MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.([1.0], [x]), 0.0), MOI.GreaterThan(1.0))
    c2 = MOI.add_constraint(model, MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.([1.0], [x]), 0.0), MOI.LessThan(2.0))

    # Getting the results before the conflict refiner has been called must return an error.
    @test MOI.get(model, Gurobi.ConflictStatus()) == -1
    @test MOI.get(model, MOI.ConflictStatus()) == MOI.COMPUTE_CONFLICT_NOT_CALLED
    @test_throws ErrorException MOI.get(model, MOI.ConstraintConflictStatus(), c1)

    # TODO(odow): bypass Gurobi's IIS checker when underlying model is
    # feasible.
    # @test_throws Gurobi.GurobiError MOI.compute_conflict!(model)
    # @test MOI.get(model, Gurobi.ConflictStatus()) == Gurobi.IIS_NOT_INFEASIBLE
    # @test MOI.get(model, MOI.ConflictStatus()) == MOI.NO_CONFLICT_EXISTS
    # @test MOI.get(model, MOI.ConstraintConflictStatus(), c1) == MOI.NOT_IN_CONFLICT
    # @test MOI.get(model, MOI.ConstraintConflictStatus(), c2) == MOI.NOT_IN_CONFLICT
end

function test_RawParameter()
    model = Gurobi.Optimizer(GRB_ENV)
    @test MOI.get(model, MOI.RawParameter("OutputFlag")) == 1
    MOI.set(model, MOI.RawParameter("OutputFlag"), 0)
    @test MOI.get(model, MOI.RawParameter("OutputFlag")) == 0
end

function test_QCPDual_1()
    model = Gurobi.Optimizer(GRB_ENV)
    MOI.set(model, MOI.Silent(), true)
    MOI.set(model, MOI.RawParameter("QCPDual"), 1)
    MOI.Utilities.loadfromstring!(model, """
    variables: x, y, z
    minobjective: 1.0 * x + 1.0 * y + 1.0 * z
    c1: x + y == 2.0
    c2: x + y + z >= 0.0
    c3: 1.0 * x * x + -1.0 * y * y + -1.0 * z * z >= 0.0
    c4: x >= 0.0
    c5: y >= 0.0
    c6: z >= 0.0
    """)
    MOI.optimize!(model)
    @test MOI.get(model, MOI.TerminationStatus()) == MOI.OPTIMAL
    @test MOI.get(model, MOI.PrimalStatus()) == MOI.FEASIBLE_POINT
    @test MOI.get(model, MOI.DualStatus()) == MOI.FEASIBLE_POINT
    c1 = MOI.get(model, MOI.ConstraintIndex, "c1")
    c2 = MOI.get(model, MOI.ConstraintIndex, "c2")
    c3 = MOI.get(model, MOI.ConstraintIndex, "c3")
    @test MOI.get(model, MOI.ConstraintDual(), c1) ≈ 1.0 atol=1e-6
    @test MOI.get(model, MOI.ConstraintDual(), c2) ≈ 0.0 atol=1e-6
    @test MOI.get(model, MOI.ConstraintDual(), c3) ≈ 0.0 atol=1e-6
end

function test_QCPDual_default()
    model = Gurobi.Optimizer(GRB_ENV)
    MOI.set(model, MOI.Silent(), true)
    MOI.Utilities.loadfromstring!(model, """
    variables: x, y, z
    minobjective: 1.0 * x + 1.0 * y + 1.0 * z
    c1: x + y == 2.0
    c2: x + y + z >= 0.0
    c3: 1.0 * x * x + -1.0 * y * y + -1.0 * z * z >= 0.0
    c4: x >= 0.0
    c5: y >= 0.0
    c6: z >= 0.0
    """)
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
end

function test_Add_and_delete_constraints()
    model = Gurobi.Optimizer(GRB_ENV)
    x = MOI.add_variables(model, 2)
    cs = MOI.add_constraints(
        model,
        [MOI.ScalarAffineFunction([MOI.ScalarAffineTerm(1.0, x[i])], 0.0) for i in 1:2],
        MOI.EqualTo.([0.0, 0.0])
    )
    @test MOI.get(model, MOI.NumberOfConstraints{
        MOI.ScalarAffineFunction{Float64}, MOI.EqualTo{Float64}
    }()) == 2
    MOI.delete(model, cs)
    @test iszero(MOI.get(model, MOI.NumberOfConstraints{
        MOI.ScalarAffineFunction{Float64}, MOI.EqualTo{Float64}
    }()))
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
end


function test_extra_name_Variables()
    model = Gurobi.Optimizer(GRB_ENV)
    x = MOI.add_variables(model, 3)
    MOI.set(model, MOI.VariableName(), x[1], "x1")
    @test MOI.get(model, MOI.VariableIndex, "x1") == x[1]
    MOI.set(model, MOI.VariableName(), x[1], "x2")
    @test MOI.get(model, MOI.VariableIndex, "x1") === nothing
    @test MOI.get(model, MOI.VariableIndex, "x2") == x[1]
    MOI.set(model, MOI.VariableName(), x[2], "x1")
    @test MOI.get(model, MOI.VariableIndex, "x1") == x[2]
    MOI.set(model, MOI.VariableName(), x[3], "xα")
    @test MOI.get(model, MOI.VariableIndex, "xα") == x[3]
    MOI.set(model, MOI.VariableName(), x[1], "x1")
    @test_throws ErrorException MOI.get(model, MOI.VariableIndex, "x1")
end

function test_extra_name_Variable_bounds()
    model = Gurobi.Optimizer(GRB_ENV)
    x = MOI.add_variable(model)
    c1 = MOI.add_constraint(model, MOI.SingleVariable(x), MOI.GreaterThan(0.0))
    c2 = MOI.add_constraint(model, MOI.SingleVariable(x), MOI.LessThan(1.0))
    MOI.set(model, MOI.ConstraintName(), c1, "c1")
    @test MOI.get(model, MOI.ConstraintIndex, "c1") == c1
    MOI.set(model, MOI.ConstraintName(), c1, "c2")
    @test MOI.get(model, MOI.ConstraintIndex, "c1") === nothing
    @test MOI.get(model, MOI.ConstraintIndex, "c2") == c1
    MOI.set(model, MOI.ConstraintName(), c2, "c1")
    @test MOI.get(model, MOI.ConstraintIndex, "c1") == c2
    MOI.set(model, MOI.ConstraintName(), c1, "c1")
    @test_throws ErrorException MOI.get(model, MOI.ConstraintIndex, "c1")
end

function test_extra_name_Affine_constraints()
    model = Gurobi.Optimizer(GRB_ENV)
    x = MOI.add_variable(model)
    f = MOI.ScalarAffineFunction([MOI.ScalarAffineTerm(1.0, x)], 0.0)
    c1 = MOI.add_constraint(model, f, MOI.GreaterThan(0.0))
    c2 = MOI.add_constraint(model, f, MOI.LessThan(1.0))
    MOI.set(model, MOI.ConstraintName(), c1, "c1")
    @test MOI.get(model, MOI.ConstraintIndex, "c1") == c1
    MOI.set(model, MOI.ConstraintName(), c1, "c2")
    @test MOI.get(model, MOI.ConstraintIndex, "c1") === nothing
    @test MOI.get(model, MOI.ConstraintIndex, "c2") == c1
    MOI.set(model, MOI.ConstraintName(), c2, "c1")
    @test MOI.get(model, MOI.ConstraintIndex, "c1") == c2
    MOI.set(model, MOI.ConstraintName(), c1, "c1")
    @test_throws ErrorException MOI.get(model, MOI.ConstraintIndex, "c1")
end

function test_ConstraintAttribute()
    model = Gurobi.Optimizer(GRB_ENV)
    MOI.Utilities.loadfromstring!(model, """
variables: x
minobjective: x
c1: x >= 0.0
c2: 2x >= 1.0
c3: x in Integer()
""")
    c2 = MOI.get(model, MOI.ConstraintIndex, "c2")
    c3 = MOI.get(model, MOI.ConstraintIndex, "c3")
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
    # Things that should fail follow.
    # Non-linear constraints are not supported.
    @test_throws(
        MOI.SetAttributeNotAllowed(Gurobi.ConstraintAttribute("Lazy")),
        MOI.set(model, Gurobi.ConstraintAttribute("Lazy"), c3, 1)
    )
    # Getting/setting a non-existing attribute.
    attr = Gurobi.ConstraintAttribute("Non-existing")
    err = ErrorException("Gurobi Error 10004: Unknown attribute 'Non-existing'")
    @test_throws err MOI.set(model, attr, c2, 1)
    @test_throws err MOI.get(model, attr, c2)
    # Setting an attribute to a value of the wrong type.
    @test_throws(
        ArgumentError("Attribute Lazy requires Int64 arguments. Provided argument was of type Float64."),
        MOI.set(model, Gurobi.ConstraintAttribute("Lazy"), c2, 1.5)
    )
end

function test_VariableAttribute()
    model = Gurobi.Optimizer(GRB_ENV)
    MOI.Utilities.loadfromstring!(model, """
variables: x
minobjective: x
c1: x >= 0.0
c2: 2x >= 1.0
c3: x in Integer()
""")
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
        ArgumentError("Attribute BranchPriority requires Int64 arguments. Provided argument was of type Float64."),
        MOI.set(model, Gurobi.VariableAttribute("BranchPriority"), x, 1.5)
    )
end

function test_ModelAttribute()
    model = Gurobi.Optimizer(GRB_ENV)
    MOI.Utilities.loadfromstring!(model, """
variables: x
minobjective: x
c1: x >= 0.0
c2: 2x >= 1.0
c3: x in Integer()
""")
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
        ArgumentError("Attribute NumStart requires $(Int) arguments. Provided argument was of type Float64."),
        MOI.set(model, Gurobi.ModelAttribute("NumStart"), 4.0)
    )
end

function test_soc_no_initial_bound()
    model = Gurobi.Optimizer(GRB_ENV)
    MOI.set(model, MOI.Silent(), true)
    t = MOI.add_variable(model)
    x = MOI.add_variables(model, 2)
    MOI.add_constraints(
        model, MOI.SingleVariable.(x), MOI.GreaterThan.(3.0:4.0)
    )
    c_soc = MOI.add_constraint(
        model, MOI.VectorOfVariables([t; x]), MOI.SecondOrderCone(3)
    )
    MOI.set(model, MOI.ObjectiveFunction{MOI.SingleVariable}(), MOI.SingleVariable(t))
    MOI.set(model, MOI.ObjectiveSense(), MOI.MIN_SENSE)
    MOI.optimize!(model)
    @test MOI.get(model, MOI.VariablePrimal(), t) == 5.0
    MOI.delete(model, c_soc)
    MOI.optimize!(model)
    @test MOI.get(model, MOI.TerminationStatus()) == MOI.DUAL_INFEASIBLE
end

function test_soc_nonnegative_initial_bound()
    model = Gurobi.Optimizer(GRB_ENV)
    MOI.set(model, MOI.Silent(), true)
    t = MOI.add_variable(model)
    x = MOI.add_variables(model, 2)
    MOI.add_constraint(model, MOI.SingleVariable(t), MOI.GreaterThan(1.0))
    MOI.add_constraints(
        model, MOI.SingleVariable.(x), MOI.GreaterThan.(3.0:4.0)
    )
    c_soc = MOI.add_constraint(
        model, MOI.VectorOfVariables([t; x]), MOI.SecondOrderCone(3)
    )
    MOI.set(model, MOI.ObjectiveFunction{MOI.SingleVariable}(), MOI.SingleVariable(t))
    MOI.set(model, MOI.ObjectiveSense(), MOI.MIN_SENSE)
    MOI.optimize!(model)
    @test MOI.get(model, MOI.VariablePrimal(), t) == 5.0
    MOI.delete(model, c_soc)
    MOI.optimize!(model)
    @test MOI.get(model, MOI.VariablePrimal(), t) == 1.0
end

function test_soc_negative_initial_bound()
    model = Gurobi.Optimizer(GRB_ENV)
    MOI.set(model, MOI.Silent(), true)
    t = MOI.add_variable(model)
    x = MOI.add_variables(model, 2)
    MOI.add_constraint(model, MOI.SingleVariable(t), MOI.GreaterThan(-1.0))
    MOI.add_constraints(
        model, MOI.SingleVariable.(x), MOI.GreaterThan.(3.0:4.0)
    )
    c_soc = MOI.add_constraint(
        model, MOI.VectorOfVariables([t; x]), MOI.SecondOrderCone(3)
    )
    MOI.optimize!(model)
    MOI.set(model, MOI.ObjectiveFunction{MOI.SingleVariable}(), MOI.SingleVariable(t))
    MOI.set(model, MOI.ObjectiveSense(), MOI.MIN_SENSE)
    @test MOI.get(model, MOI.VariablePrimal(), t) == 5.0
    MOI.delete(model, c_soc)
    MOI.optimize!(model)
    @test MOI.get(model, MOI.VariablePrimal(), t) == -1.0
end

function test_soc_nonnegative_post_bound()
    model = Gurobi.Optimizer(GRB_ENV)
    MOI.set(model, MOI.Silent(), true)
    t = MOI.add_variable(model)
    x = MOI.add_variables(model, 2)
    MOI.add_constraints(
        model, MOI.SingleVariable.(x), MOI.GreaterThan.(3.0:4.0)
    )
    MOI.add_constraint(
        model, MOI.VectorOfVariables([t; x]), MOI.SecondOrderCone(3)
    )
    c_lb = MOI.add_constraint(model, MOI.SingleVariable(t), MOI.GreaterThan(6.0))
    MOI.set(model, MOI.ObjectiveFunction{MOI.SingleVariable}(), MOI.SingleVariable(t))
    MOI.set(model, MOI.ObjectiveSense(), MOI.MIN_SENSE)
    MOI.optimize!(model)
    @test MOI.get(model, MOI.VariablePrimal(), t) == 6.0
    MOI.delete(model, c_lb)
    MOI.optimize!(model)
    @test MOI.get(model, MOI.VariablePrimal(), t) == 5.0
end

function test_soc_negative_post_bound()
    model = Gurobi.Optimizer(GRB_ENV)
    MOI.set(model, MOI.Silent(), true)
    t = MOI.add_variable(model)
    x = MOI.add_variables(model, 2)
    MOI.add_constraints(
        model, MOI.SingleVariable.(x), MOI.GreaterThan.(3.0:4.0)
    )
    c_soc = MOI.add_constraint(
        model, MOI.VectorOfVariables([t; x]), MOI.SecondOrderCone(3)
    )
    MOI.add_constraint(model, MOI.SingleVariable(t), MOI.GreaterThan(-6.0))
    MOI.set(model, MOI.ObjectiveFunction{MOI.SingleVariable}(), MOI.SingleVariable(t))
    MOI.set(model, MOI.ObjectiveSense(), MOI.MIN_SENSE)
    MOI.optimize!(model)
    @test MOI.get(model, MOI.VariablePrimal(), t) == 5.0
    MOI.delete(model, c_soc)
    MOI.optimize!(model)
    @test MOI.get(model, MOI.VariablePrimal(), t) == -6.0
end

function test_soc_negative_post_bound_ii()
    model = Gurobi.Optimizer(GRB_ENV)
    MOI.set(model, MOI.Silent(), true)
    t = MOI.add_variable(model)
    x = MOI.add_variables(model, 2)
    MOI.add_constraints(
        model, MOI.SingleVariable.(x), MOI.GreaterThan.(3.0:4.0)
    )
    c_soc = MOI.add_constraint(
        model, MOI.VectorOfVariables([t; x]), MOI.SecondOrderCone(3)
    )
    c_lb = MOI.add_constraint(model, MOI.SingleVariable(t), MOI.GreaterThan(-6.0))
    MOI.set(model, MOI.ObjectiveFunction{MOI.SingleVariable}(), MOI.SingleVariable(t))
    MOI.set(model, MOI.ObjectiveSense(), MOI.MIN_SENSE)
    MOI.optimize!(model)
    @test MOI.get(model, MOI.VariablePrimal(), t) == 5.0
    MOI.delete(model, c_lb)
    MOI.optimize!(model)
    @test MOI.get(model, MOI.VariablePrimal(), t) == 5.0
    MOI.delete(model, c_soc)
    MOI.optimize!(model)
    @test MOI.get(model, MOI.TerminationStatus()) == MOI.DUAL_INFEASIBLE
end

function test_soc_negative_post_bound_iii()
    # This test was added because add_constraint and add_constraints had
    # different implementations and add_constraints failed where
    # add_constraint succeeded.
    model = Gurobi.Optimizer(GRB_ENV)
    MOI.set(model, MOI.Silent(), true)
    t = MOI.add_variable(model)
    x = MOI.add_variables(model, 2)
    c_soc = MOI.add_constraint(
        model, MOI.VectorOfVariables([t; x]), MOI.SecondOrderCone(3)
    )
    c_lbs = MOI.add_constraints(
        model, MOI.SingleVariable.([t; x]), MOI.GreaterThan.([-6.0, 3.0, 4.0])
    )
    c_lb = first(c_lbs)
    MOI.set(model, MOI.ObjectiveFunction{MOI.SingleVariable}(), MOI.SingleVariable(t))
    MOI.set(model, MOI.ObjectiveSense(), MOI.MIN_SENSE)
    MOI.optimize!(model)
    @test MOI.get(model, MOI.VariablePrimal(), t) == 5.0
    MOI.delete(model, c_lb)
    MOI.optimize!(model)
    @test MOI.get(model, MOI.VariablePrimal(), t) == 5.0
    MOI.delete(model, c_soc)
    MOI.optimize!(model)
    @test MOI.get(model, MOI.TerminationStatus()) == MOI.DUAL_INFEASIBLE
end

function test_2_soc()
    model = Gurobi.Optimizer(GRB_ENV)
    MOI.set(model, MOI.Silent(), true)
    t = MOI.add_variable(model)
    x = MOI.add_variables(model, 2)
    MOI.add_constraint(model, MOI.SingleVariable(t), MOI.GreaterThan(-1.0))
    MOI.add_constraints(
        model, MOI.SingleVariable.(x), MOI.GreaterThan.([4.0, 3.0])
    )
    c_soc_1 = MOI.add_constraint(
        model, MOI.VectorOfVariables([t, x[1]]), MOI.SecondOrderCone(2)
    )
    c_soc_2 = MOI.add_constraint(
        model, MOI.VectorOfVariables([t, x[2]]), MOI.SecondOrderCone(2)
    )
    MOI.optimize!(model)
    MOI.set(model, MOI.ObjectiveFunction{MOI.SingleVariable}(), MOI.SingleVariable(t))
    MOI.set(model, MOI.ObjectiveSense(), MOI.MIN_SENSE)
    @test MOI.get(model, MOI.VariablePrimal(), t) == 4.0
    MOI.delete(model, c_soc_1)
    MOI.optimize!(model)
    @test MOI.get(model, MOI.VariablePrimal(), t) == 3.0
    MOI.delete(model, c_soc_2)
    MOI.optimize!(model)
    @test MOI.get(model, MOI.VariablePrimal(), t) == -1.0
end

function test_Duplicate_names_Variables()
    model = Gurobi.Optimizer(GRB_ENV)
    (x, y, z) = MOI.add_variables(model, 3)
    MOI.set(model, MOI.VariableName(), x, "x")
    MOI.set(model, MOI.VariableName(), y, "x")
    MOI.set(model, MOI.VariableName(), z, "z")
    @test MOI.get(model, MOI.VariableIndex, "z") == z
    @test_throws ErrorException MOI.get(model, MOI.VariableIndex, "x")
    MOI.set(model, MOI.VariableName(), y, "y")
    @test MOI.get(model, MOI.VariableIndex, "x") == x
    @test MOI.get(model, MOI.VariableIndex, "y") == y
    MOI.set(model, MOI.VariableName(), z, "x")
    @test_throws ErrorException MOI.get(model, MOI.VariableIndex, "x")
    MOI.delete(model, x)
    @test MOI.get(model, MOI.VariableIndex, "x") == z
end

function test_Duplicate_names_SingleVariable()
    model = Gurobi.Optimizer(GRB_ENV)
    x = MOI.add_variables(model, 3)
    c = MOI.add_constraints(model, MOI.SingleVariable.(x), MOI.GreaterThan(0.0))
    MOI.set(model, MOI.ConstraintName(), c[1], "x")
    MOI.set(model, MOI.ConstraintName(), c[2], "x")
    MOI.set(model, MOI.ConstraintName(), c[3], "z")
    @test MOI.get(model, MOI.ConstraintIndex, "z") == c[3]
    @test_throws ErrorException MOI.get(model, MOI.ConstraintIndex, "x")
    MOI.set(model, MOI.ConstraintName(), c[2], "y")
    @test MOI.get(model, MOI.ConstraintIndex, "x") == c[1]
    @test MOI.get(model, MOI.ConstraintIndex, "y") == c[2]
    MOI.set(model, MOI.ConstraintName(), c[3], "x")
    @test_throws ErrorException MOI.get(model, MOI.ConstraintIndex, "x")
    MOI.delete(model, c[1])
    @test MOI.get(model, MOI.ConstraintIndex, "x") == c[3]
    MOI.set(model, MOI.ConstraintName(), c[2], "x")
    @test_throws ErrorException MOI.get(model, MOI.ConstraintIndex, "x")
    MOI.delete(model, x[3])
    @test MOI.get(model, MOI.ConstraintIndex, "x") == c[2]
end

function test_Duplicate_names_ScalarAffineFunction()
    model = Gurobi.Optimizer(GRB_ENV)
    x = MOI.add_variables(model, 3)
    fs = [
        MOI.ScalarAffineFunction([MOI.ScalarAffineTerm(1.0, xi)], 0.0)
        for xi in x
    ]
    c = MOI.add_constraints(model, fs, MOI.GreaterThan(0.0))
    MOI.set(model, MOI.ConstraintName(), c[1], "x")
    MOI.set(model, MOI.ConstraintName(), c[2], "x")
    MOI.set(model, MOI.ConstraintName(), c[3], "z")
    @test MOI.get(model, MOI.ConstraintIndex, "z") == c[3]
    @test_throws ErrorException MOI.get(model, MOI.ConstraintIndex, "x")
    MOI.set(model, MOI.ConstraintName(), c[2], "y")
    @test MOI.get(model, MOI.ConstraintIndex, "x") == c[1]
    @test MOI.get(model, MOI.ConstraintIndex, "y") == c[2]
    MOI.set(model, MOI.ConstraintName(), c[3], "x")
    @test_throws ErrorException MOI.get(model, MOI.ConstraintIndex, "x")
    MOI.delete(model, c[1])
    @test MOI.get(model, MOI.ConstraintIndex, "x") == c[3]
end

function test_Duals_with_equal_bounds_250()
    model = Gurobi.Optimizer(GRB_ENV)
    MOI.set(model, MOI.Silent(), true)
    x = MOI.add_variable(model)
    xl = MOI.add_constraint(model, MOI.SingleVariable(x), MOI.GreaterThan(1.0))
    xu = MOI.add_constraint(model, MOI.SingleVariable(x), MOI.LessThan(1.0))
    MOI.set(model, MOI.ObjectiveSense(), MOI.MIN_SENSE)
    MOI.set(model, MOI.ObjectiveFunction{MOI.SingleVariable}(), MOI.SingleVariable(x))
    MOI.optimize!(model)
    @test MOI.get(model, MOI.ConstraintDual(), xl) == 1.0
    @test MOI.get(model, MOI.ConstraintDual(), xu) == 0.0
end

function test_Objective_functions()
    model = Gurobi.Optimizer(GRB_ENV)
    x = MOI.add_variable(model)
    @test MOI.get(model, MOI.ObjectiveSense()) == MOI.FEASIBILITY_SENSE
    @test MOI.get(model, MOI.ListOfModelAttributesSet()) == Any[MOI.ObjectiveSense()]
    MOI.set(model, MOI.ObjectiveSense(), MOI.MIN_SENSE)
    MOI.set(model, MOI.ObjectiveFunction{MOI.SingleVariable}(), MOI.SingleVariable(x))
    @test MOI.get(model, MOI.ListOfModelAttributesSet()) ==
        Any[MOI.ObjectiveSense(), MOI.ObjectiveFunction{MOI.SingleVariable}()]
    MOI.set(model, MOI.ObjectiveSense(), MOI.FEASIBILITY_SENSE)
    @test MOI.get(model, MOI.ListOfModelAttributesSet()) == Any[MOI.ObjectiveSense()]
end

function test_FEASIBILITY_SENSE_zeros_objective()
    model = Gurobi.Optimizer(GRB_ENV)
    MOI.set(model, MOI.Silent(), true)
    x = MOI.add_variable(model)
    MOI.add_constraint(model, MOI.SingleVariable(x), MOI.GreaterThan(1.0))
    MOI.set(model, MOI.ObjectiveSense(), MOI.MIN_SENSE)
    MOI.set(model, MOI.ObjectiveFunction{MOI.SingleVariable}(), MOI.SingleVariable(x))
    MOI.optimize!(model)
    @test MOI.get(model, MOI.ObjectiveValue()) == 1.0
    MOI.set(model, MOI.ObjectiveSense(), MOI.FEASIBILITY_SENSE)
    MOI.optimize!(model)
    @test MOI.get(model, MOI.ObjectiveValue()) == 0.0
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
    MOI.set(model, Gurobi.CallbackFunction(), (cb_data, cb_where) -> begin
        GRBterminate(model)
    end)
    MOI.optimize!(model)
    @test MOI.get(model, MOI.TerminationStatus()) == MOI.INTERRUPTED
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
    MOI.add_constraint(model, MOI.SingleVariable(x), MOI.Integer())
    MOI.set(model, MOI.ObjectiveSense(), MOI.MIN_SENSE)
    MOI.set(
        model,
        MOI.ObjectiveFunction{MOI.SingleVariable}(),
        MOI.SingleVariable(x),
    )
    MOI.set(model, MOI.RawParameter("LazyConstraints"), 1)
    interrupt_thrown = false
    i = 0.0
    MOI.set(model, Gurobi.CallbackFunction(), (cb_data, cb_where) -> begin
        if cb_where != Gurobi.GRB_CB_MIPSOL
            return
        end
        MOI.submit(
            model,
            MOI.LazyConstraint(cb_data),
            MOI.ScalarAffineFunction{Float64}(
                [MOI.ScalarAffineTerm(1.0, x)], 0.0
            ),
            MOI.GreaterThan{Float64}(i)
        )
        i += 1
        if !interrupt_thrown
            interrupt_thrown = true
            throw(InterruptException())
        end
    end)
    MOI.optimize!(model)
    @test MOI.get(model, MOI.TerminationStatus()) == MOI.INTERRUPTED
end

function test_indicator_name()
    MOI.empty!(OPTIMIZER)
    x = MOI.add_variables(model, 2)
    MOI.add_constraint(model, MOI.SingleVariable(x[1]), MOI.ZeroOne())
    f = MOI.VectorAffineFunction(
        [
            MOI.VectorAffineTerm(1, MOI.ScalarAffineTerm(1.0, x[1])),
            MOI.VectorAffineTerm(2, MOI.ScalarAffineTerm(2.0, x[2])),
        ],
        [0.0, 0.0],
    )
    s = MOI.IndicatorSet{MOI.ACTIVATE_ON_ONE}(MOI.GreaterThan(1.0))
    c = MOI.add_constraint(model, f, s)
    MOI.set(model, MOI.ConstraintName(), c, "my_indicator")
    @test MOI.get(model, MOI.ConstraintName(), c) == "my_indicator"
end

function test_indicator_on_one()
    MOI.empty!(OPTIMIZER)
    x = MOI.add_variables(model, 2)
    MOI.add_constraint(model, MOI.SingleVariable(x[1]), MOI.ZeroOne())
    f = MOI.VectorAffineFunction(
        [
            MOI.VectorAffineTerm(1, MOI.ScalarAffineTerm(1.0, x[1])),
            MOI.VectorAffineTerm(2, MOI.ScalarAffineTerm(2.0, x[2])),
        ],
        [0.0, 0.0],
    )
    s = MOI.IndicatorSet{MOI.ACTIVATE_ON_ONE}(MOI.GreaterThan(1.0))
    c = MOI.add_constraint(model, f, s)
    @test MOI.get(model, MOI.ConstraintSet(), c) == s
    @test isapprox(MOI.get(model, MOI.ConstraintFunction(), c), f)
end

function test_indicator_on_zero()
    MOI.empty!(OPTIMIZER)
    x = MOI.add_variables(model, 2)
    MOI.add_constraint(model, MOI.SingleVariable(x[1]), MOI.ZeroOne())
    f = MOI.VectorAffineFunction(
        [
            MOI.VectorAffineTerm(1, MOI.ScalarAffineTerm(1.0, x[1])),
            MOI.VectorAffineTerm(2, MOI.ScalarAffineTerm(2.0, x[2])),
        ],
        [0.0, 0.0],
    )
    s = MOI.IndicatorSet{MOI.ACTIVATE_ON_ZERO}(MOI.GreaterThan(1.0))
    c = MOI.add_constraint(model, f, s)
    @test MOI.get(model, MOI.ConstraintSet(), c) == s
    @test isapprox(MOI.get(model, MOI.ConstraintFunction(), c), f)
end

function test_indicator_nonconstant_x()
    MOI.empty!(OPTIMIZER)
    x = MOI.add_variables(model, 2)
    MOI.add_constraint(model, MOI.SingleVariable(x[1]), MOI.ZeroOne())
    f = MOI.VectorAffineFunction(
        [
            MOI.VectorAffineTerm(1, MOI.ScalarAffineTerm(-1.0, x[1])),
            MOI.VectorAffineTerm(2, MOI.ScalarAffineTerm(2.0, x[2])),
        ],
        [0.0, 0.0],
    )
    s = MOI.IndicatorSet{MOI.ACTIVATE_ON_ONE}(MOI.GreaterThan(1.0))
    @test_throws ErrorException MOI.add_constraint(model, f, s)
end

function test_indicator_too_many_indicators()
    MOI.empty!(OPTIMIZER)
    x = MOI.add_variables(model, 2)
    MOI.add_constraint(model, MOI.SingleVariable(x[1]), MOI.ZeroOne())
    f = MOI.VectorAffineFunction(
        [
            MOI.VectorAffineTerm(1, MOI.ScalarAffineTerm(-1.0, x[1])),
            MOI.VectorAffineTerm(1, MOI.ScalarAffineTerm(2.0, x[2])),
        ],
        [0.0, 0.0],
    )
    s = MOI.IndicatorSet{MOI.ACTIVATE_ON_ONE}(MOI.GreaterThan(1.0))
    @test_throws ErrorException MOI.add_constraint(model, f, s)
end

function test_indicator_nonconstant()
    MOI.empty!(OPTIMIZER)
    x = MOI.add_variables(model, 2)
    MOI.add_constraint(model, MOI.SingleVariable(x[1]), MOI.ZeroOne())
    f = MOI.VectorAffineFunction(
        [
            MOI.VectorAffineTerm(1, MOI.ScalarAffineTerm(1.0, x[1])),
            MOI.VectorAffineTerm(2, MOI.ScalarAffineTerm(2.0, x[2])),
        ],
        [1.0, 0.0],
    )
    s = MOI.IndicatorSet{MOI.ACTIVATE_ON_ONE}(MOI.GreaterThan(1.0))
    @test_throws ErrorException MOI.add_constraint(model, f, s)
end

end

runtests(TestMOIWrapper)
