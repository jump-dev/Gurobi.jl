const MOI  = Gurobi.MOI
const MOIT = MOI.Test

const GUROBI_ENV = Gurobi.Env()
const OPTIMIZER = MOI.Bridges.full_bridge_optimizer(
    Gurobi.Optimizer(GUROBI_ENV, OutputFlag=0), Float64
)

const CONFIG = MOIT.TestConfig()

@testset "Unit Tests" begin
    MOIT.basic_constraint_tests(OPTIMIZER, CONFIG)
    MOIT.unittest(OPTIMIZER, MOIT.TestConfig(atol=1e-6))
    MOIT.modificationtest(OPTIMIZER, CONFIG)
end

@testset "Linear tests" begin
    @testset "Default Solver"  begin
        MOIT.contlineartest(OPTIMIZER, MOIT.TestConfig(basis = true), [
            # This requires an infeasiblity certificate for a variable bound.
            "linear12"
        ])
    end
    @testset "No certificate" begin
        MOIT.linear12test(OPTIMIZER, MOIT.TestConfig(infeas_certificates=false))
    end
end

@testset "Quadratic tests" begin
    MOIT.contquadratictest(OPTIMIZER, MOIT.TestConfig(atol=1e-3, rtol=1e-3), [
        "ncqcp"  # Gurobi doesn't support non-convex problems.
    ])
end

@testset "Linear Conic tests" begin
    MOIT.lintest(OPTIMIZER, CONFIG)
end

@testset "Integer Linear tests" begin
    MOIT.intlineartest(OPTIMIZER, CONFIG, [
        # Indicator sets not supported.
        "indicator1", "indicator2", "indicator3"
    ])
end

@testset "ModelLike tests" begin

    @test MOI.get(OPTIMIZER, MOI.SolverName()) == "Gurobi"

    @testset "default_objective_test" begin
        MOIT.default_objective_test(OPTIMIZER)
    end

    @testset "default_status_test" begin
        MOIT.default_status_test(OPTIMIZER)
    end

    @testset "nametest" begin
        MOIT.nametest(OPTIMIZER)
    end

    @testset "validtest" begin
        MOIT.validtest(OPTIMIZER)
    end

    @testset "emptytest" begin
        MOIT.emptytest(OPTIMIZER)
    end

    @testset "orderedindicestest" begin
        MOIT.orderedindicestest(OPTIMIZER)
    end

    @testset "copytest" begin
        MOIT.copytest(
            OPTIMIZER,
            MOI.Bridges.full_bridge_optimizer(Gurobi.Optimizer(GUROBI_ENV), Float64)
        )
    end

    @testset "scalar_function_constant_not_zero" begin
        MOIT.scalar_function_constant_not_zero(OPTIMIZER)
    end

    @testset "start_values_test" begin
        model = Gurobi.Optimizer(GUROBI_ENV, OutputFlag = 0)
        x = MOI.add_variables(model, 2)
        @test MOI.supports(model, MOI.VariablePrimalStart(), MOI.VariableIndex)
        MOI.set(model, MOI.VariablePrimalStart(), x[1], 1.0)
        MOI.set(model, MOI.VariablePrimalStart(), x[2], nothing)
        @test MOI.get(model, MOI.VariablePrimalStart(), x[1]) == 1.0
        @test MOI.get(model, MOI.VariablePrimalStart(), x[2]) === nothing
        MOI.optimize!(model)
        @test MOI.get(model, MOI.ObjectiveValue()) == 0.0
        # We don't support ConstraintDualStart or ConstraintPrimalStart yet.
        # @test_broken MOIT.start_values_test(Gurobi.Optimizer(GUROBI_ENV), OPTIMIZER)
    end

    @testset "supports_constrainttest" begin
        # supports_constrainttest needs VectorOfVariables-in-Zeros,
        # MOIT.supports_constrainttest(Gurobi.Optimizer(GUROBI_ENV), Float64, Float32)
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

    @testset "set_lower_bound_twice" begin
        MOIT.set_lower_bound_twice(OPTIMIZER, Float64)
    end

    @testset "set_upper_bound_twice" begin
        MOIT.set_upper_bound_twice(OPTIMIZER, Float64)
    end
end

@testset "Gurobi Callback" begin
    @testset "Generic callback" begin
        m = Gurobi.Optimizer(GUROBI_ENV, OutputFlag=0)
        x = MOI.add_variable(m)
        MOI.add_constraint(m, MOI.SingleVariable(x), MOI.GreaterThan(1.0))
        MOI.set(m, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(),
            MOI.ScalarAffineFunction{Float64}(
                [MOI.ScalarAffineTerm{Float64}(1.0, x)],
                0.0
            )
        )

        cb_calls = Int32[]
        function callback_function(cb_data::Gurobi.CallbackData, cb_where::Int32)
            push!(cb_calls, cb_where)
            nothing
        end

        MOI.set(m, Gurobi.CallbackFunction(), callback_function)
        MOI.optimize!(m)

        @test length(cb_calls) > 0
        @test Gurobi.CB_MESSAGE in cb_calls
        @test Gurobi.CB_PRESOLVE in cb_calls
        @test !(Gurobi.CB_MIPSOL in cb_calls)
    end

    @testset "Lazy cut" begin
        m = Gurobi.Optimizer(GUROBI_ENV, OutputFlag=0, Cuts=0, Presolve=0, Heuristics=0, LazyConstraints=1)
        MOI.Utilities.loadfromstring!(m,"""
            variables: x, y
            maxobjective: y
            c1: x in Integer()
            c2: y in Integer()
            c3: x in Interval(0.0, 2.0)
            c4: y in Interval(0.0, 2.0)
        """)
        x = MOI.get(m, MOI.VariableIndex, "x")
        y = MOI.get(m, MOI.VariableIndex, "y")

        # We now define our callback function that takes two arguments:
        #   (1) the callback handle; and
        #   (2) the location from where the callback was called.
        # Note that we can access m, x, and y because this function is defined
        # inside the same scope
        cb_calls = Int32[]
        function callback_function(cb_data::Gurobi.CallbackData, cb_where::Int32)
            push!(cb_calls, cb_where)
            if cb_where == Gurobi.CB_MIPSOL
                Gurobi.load_callback_variable_primal(m, cb_data, cb_where)
                x_val = MOI.get(m, Gurobi.CallbackVariablePrimal(), x)
                y_val = MOI.get(m, Gurobi.CallbackVariablePrimal(), y)
                # We have two constraints, one cutting off the top
                # left corner and one cutting off the top right corner, e.g.
                # (0,2) +---+---+ (2,2)
                #       |xx/ \xx|
                #       |x/   \x|
                #       |/     \|
                # (0,1) +       + (2,1)
                #       |       |
                # (0,0) +---+---+ (2,0)
                TOL = 1e-6  # Allow for some impreciseness in the solution
                if y_val - x_val > 1 + TOL
                    Gurobi.cblazy!(cb_data, m,
                        MOI.ScalarAffineFunction{Float64}(
                            MOI.ScalarAffineTerm.([-1.0, 1.0], [x, y]),
                            0.0
                        ),
                        MOI.LessThan{Float64}(1.0)
                    )
                elseif y_val + x_val > 3 + TOL
                    Gurobi.cblazy!(cb_data, m,
                        MOI.ScalarAffineFunction{Float64}(
                            MOI.ScalarAffineTerm.([1.0, 1.0], [x, y]),
                            0.0
                        ),
                        MOI.LessThan{Float64}(3.0)
                    )
                end
            end
        end

        MOI.set(m, Gurobi.CallbackFunction(), callback_function)
        MOI.optimize!(m)

        @test MOI.get(m, MOI.VariablePrimal(), x) == 1
        @test MOI.get(m, MOI.VariablePrimal(), y) == 2

        @test length(cb_calls) > 0
        @test Gurobi.CB_MESSAGE in cb_calls
        @test Gurobi.CB_PRESOLVE in cb_calls
        @test Gurobi.CB_MIPSOL in cb_calls
    end
end

@testset "LQOI Issue #38" begin
    # https://github.com/JuliaOpt/LinQuadOptInterface.jl/issues/38#issuecomment-407625187
    _getinner(opt::Gurobi.Optimizer) = opt.inner
    @inferred _getinner(Gurobi.Optimizer(GUROBI_ENV))
end

@testset "User limit handling (issue #140)" begin
    # Verify that we return the correct status codes when a mixed-integer
    # problem has been solved to a *feasible* but not necessarily optimal
    # solution. To do that, we will set up an intentionally dumbed-down
    # Gurobi Gurobi.Optimizer (with all heuristics and pre-solve turned off) and
    # ask it to solve a classic knapsack problem. Setting SolutionLimit=1
    # forces the solver to return after its first feasible MIP solution,
    # which tests the right part of the code without relying on potentially
    # flaky or system-dependent time limits.
    m = Gurobi.Optimizer(GUROBI_ENV, OutputFlag=0, SolutionLimit=1,
                         Heuristics=0.0, Presolve=0)
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

@testset "Constant objective (issue #111)" begin
    m = Gurobi.Optimizer(GUROBI_ENV)
    x = MOI.add_variable(m)
    MOI.set(m, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(),
        MOI.ScalarAffineFunction(MOI.ScalarAffineTerm{Float64}[], 2.0))
    @test MOI.get(m, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}()).constant == 2.0
    @test Gurobi.get_dblattr(m.inner, "ObjCon") == 2.0

    MOI.modify(m, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(), MOI.ScalarConstantChange(3.0))
    @test MOI.get(m, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}()).constant == 3.0
    @test Gurobi.get_dblattr(m.inner, "ObjCon") == 3.0
end

@testset "Env" begin
    @testset "User-provided" begin
        env = Gurobi.Env()
        model_1 = Gurobi.Optimizer(env)
        @test model_1.inner.env === env
        model_2 = Gurobi.Optimizer(env)
        @test model_2.inner.env === env
        # Check that finalizer doesn't touch env when manually provided.
        finalize(model_1.inner)
        @test Gurobi.is_valid(env)
    end
    @testset "Automatic" begin
        model_1 = Gurobi.Optimizer()
        model_2 = Gurobi.Optimizer()
        @test model_1.inner.env !== model_2.inner.env
        # Check that env is finalized with model when not supplied manually.
        finalize(model_1.inner)
        @test !Gurobi.is_valid(model_1.inner.env)
    end
    @testset "Env when emptied" begin
        @testset "User-provided" begin
            env = Gurobi.Env()
            model = Gurobi.Optimizer(env)
            @test model.inner.env === env
            @test Gurobi.is_valid(env)
            MOI.empty!(model)
            @test model.inner.env === env
            @test Gurobi.is_valid(env)
        end
        @testset "Automatic" begin
            model = Gurobi.Optimizer()
            env = model.inner.env
            MOI.empty!(model)
            @test model.inner.env !== env
            @test Gurobi.is_valid(model.inner.env)
        end
    end
end

@testset "Conflict refiner" begin
    @testset "Variable bounds (SingleVariable and LessThan/GreaterThan)" begin
        model = Gurobi.Optimizer(GUROBI_ENV, OutputFlag=0)
        x = MOI.add_variable(model)
        c1 = MOI.add_constraint(model, MOI.SingleVariable(x), MOI.GreaterThan(2.0))
        c2 = MOI.add_constraint(model, MOI.SingleVariable(x), MOI.LessThan(1.0))

        # Getting the results before the conflict refiner has been called must return an error.
        @test MOI.get(model, Gurobi.ConflictStatus()) == MOI.OPTIMIZE_NOT_CALLED
        @test_throws ErrorException MOI.get(model, Gurobi.ConstraintConflictStatus(), c1)

        # Once it's called, no problem.
        Gurobi.compute_conflict(model)
        @test MOI.get(model, Gurobi.ConflictStatus()) == MOI.OPTIMAL
        @test MOI.get(model, Gurobi.ConstraintConflictStatus(), c1) == true
        @test MOI.get(model, Gurobi.ConstraintConflictStatus(), c2) == true
    end

    @testset "Variable bounds (ScalarAffine)" begin
        model = Gurobi.Optimizer(GUROBI_ENV, OutputFlag=0)
        x = MOI.add_variable(model)
        c1 = MOI.add_constraint(model, MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.([1.0], [x]), 0.0), MOI.GreaterThan(2.0))
        c2 = MOI.add_constraint(model, MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.([1.0], [x]), 0.0), MOI.LessThan(1.0))

        # Getting the results before the conflict refiner has been called must return an error.
        @test MOI.get(model, Gurobi.ConflictStatus()) == MOI.OPTIMIZE_NOT_CALLED
        @test_throws ErrorException MOI.get(model, Gurobi.ConstraintConflictStatus(), c1)

        # Once it's called, no problem.
        Gurobi.compute_conflict(model)
        @test MOI.get(model, Gurobi.ConflictStatus()) == MOI.OPTIMAL
        @test MOI.get(model, Gurobi.ConstraintConflictStatus(), c1) == true
        @test MOI.get(model, Gurobi.ConstraintConflictStatus(), c2) == true
    end

    @testset "Variable bounds (Invali Interval)" begin
        model = Gurobi.Optimizer(GUROBI_ENV, OutputFlag=0)
        x = MOI.add_variable(model)
        c1 = MOI.add_constraint(
            model, MOI.SingleVariable(x), MOI.Interval(1.0, 0.0)
        )
        # Getting the results before the conflict refiner has been called must return an error.
        @test MOI.get(model, Gurobi.ConflictStatus()) == MOI.OPTIMIZE_NOT_CALLED
        @test_throws ErrorException MOI.get(model, Gurobi.ConstraintConflictStatus(), c1)

        # Once it's called, no problem.
        Gurobi.compute_conflict(model)
        @test MOI.get(model, Gurobi.ConflictStatus()) == MOI.OPTIMAL
        @test MOI.get(model, Gurobi.ConstraintConflictStatus(), c1) == true
    end

    @testset "Two conflicting constraints (GreaterThan, LessThan)" begin
        model = Gurobi.Optimizer(GUROBI_ENV, OutputFlag=0)
        x = MOI.add_variable(model)
        y = MOI.add_variable(model)
        b1 = MOI.add_constraint(model, MOI.SingleVariable(x), MOI.GreaterThan(0.0))
        b2 = MOI.add_constraint(model, MOI.SingleVariable(y), MOI.GreaterThan(0.0))
        cf1 = MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.([1.0, 1.0], [x, y]), 0.0)
        c1 = MOI.add_constraint(model, cf1, MOI.LessThan(-1.0))
        cf2 = MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.([1.0, -1.0], [x, y]), 0.0)
        c2 = MOI.add_constraint(model, cf2, MOI.GreaterThan(1.0))

        # Getting the results before the conflict refiner has been called must return an error.
        @test MOI.get(model, Gurobi.ConflictStatus()) == MOI.OPTIMIZE_NOT_CALLED
        @test_throws ErrorException MOI.get(model, Gurobi.ConstraintConflictStatus(), c1)

        # Once it's called, no problem.
        Gurobi.compute_conflict(model)
        @test MOI.get(model, Gurobi.ConflictStatus()) == MOI.OPTIMAL
        @test MOI.get(model, Gurobi.ConstraintConflictStatus(), b1) == true
        @test MOI.get(model, Gurobi.ConstraintConflictStatus(), b2) == true
        @test MOI.get(model, Gurobi.ConstraintConflictStatus(), c1) == true
        @test MOI.get(model, Gurobi.ConstraintConflictStatus(), c2) == false
    end

    @testset "Two conflicting constraints (EqualTo)" begin
        model = Gurobi.Optimizer(GUROBI_ENV, OutputFlag=0)
        x = MOI.add_variable(model)
        y = MOI.add_variable(model)
        b1 = MOI.add_constraint(model, MOI.SingleVariable(x), MOI.GreaterThan(0.0))
        b2 = MOI.add_constraint(model, MOI.SingleVariable(y), MOI.GreaterThan(0.0))
        cf1 = MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.([1.0, 1.0], [x, y]), 0.0)
        c1 = MOI.add_constraint(model, cf1, MOI.EqualTo(-1.0))
        cf2 = MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.([1.0, -1.0], [x, y]), 0.0)
        c2 = MOI.add_constraint(model, cf2, MOI.GreaterThan(1.0))

        # Getting the results before the conflict refiner has been called must return an error.
        @test MOI.get(model, Gurobi.ConflictStatus()) == MOI.OPTIMIZE_NOT_CALLED
        @test_throws ErrorException MOI.get(model, Gurobi.ConstraintConflictStatus(), c1)

        # Once it's called, no problem.
        Gurobi.compute_conflict(model)
        @test MOI.get(model, Gurobi.ConflictStatus()) == MOI.OPTIMAL
        @test MOI.get(model, Gurobi.ConstraintConflictStatus(), b1) == true
        @test MOI.get(model, Gurobi.ConstraintConflictStatus(), b2) == true
        @test MOI.get(model, Gurobi.ConstraintConflictStatus(), c1) == true
        @test MOI.get(model, Gurobi.ConstraintConflictStatus(), c2) == false
    end

    @testset "Variables outside conflict" begin
        model = Gurobi.Optimizer(GUROBI_ENV, OutputFlag=0)
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
        @test MOI.get(model, Gurobi.ConflictStatus()) == MOI.OPTIMIZE_NOT_CALLED
        @test_throws ErrorException MOI.get(model, Gurobi.ConstraintConflictStatus(), c1)

        # Once it's called, no problem.
        Gurobi.compute_conflict(model)
        @test MOI.get(model, Gurobi.ConflictStatus()) == MOI.OPTIMAL
        @test MOI.get(model, Gurobi.ConstraintConflictStatus(), b1) == true
        @test MOI.get(model, Gurobi.ConstraintConflictStatus(), b2) == true
        @test MOI.get(model, Gurobi.ConstraintConflictStatus(), b3) == false
        @test MOI.get(model, Gurobi.ConstraintConflictStatus(), c1) == true
        @test MOI.get(model, Gurobi.ConstraintConflictStatus(), c2) == false
    end

    @testset "No conflict" begin
        model = Gurobi.Optimizer(GUROBI_ENV, OutputFlag=0)
        x = MOI.add_variable(model)
        c1 = MOI.add_constraint(model, MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.([1.0], [x]), 0.0), MOI.GreaterThan(1.0))
        c2 = MOI.add_constraint(model, MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.([1.0], [x]), 0.0), MOI.LessThan(2.0))

        # Getting the results before the conflict refiner has been called must return an error.
        @test MOI.get(model, Gurobi.ConflictStatus()) == MOI.OPTIMIZE_NOT_CALLED
        @test_throws ErrorException MOI.get(model, Gurobi.ConstraintConflictStatus(), c1)

        # Once it's called, no problem.
        Gurobi.compute_conflict(model)
        @test MOI.get(model, Gurobi.ConflictStatus()) == MOI.INFEASIBLE
        @test MOI.get(model, Gurobi.ConstraintConflictStatus(), c1) == false
        @test MOI.get(model, Gurobi.ConstraintConflictStatus(), c2) == false
    end
end

@testset "RawParameter" begin
    model = Gurobi.Optimizer(GUROBI_ENV)
    @test MOI.get(model, MOI.RawParameter("OutputFlag")) == 1
    MOI.set(model, MOI.RawParameter("OutputFlag"), 0)
    @test MOI.get(model, MOI.RawParameter("OutputFlag")) == 0
end

@testset "QCPDuals without needing to pass QCPDual=1" begin
    @testset "QCPDual default" begin
        model = Gurobi.Optimizer(GUROBI_ENV, OutputFlag=0)
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
    @testset "QCPDual=0" begin
        model = Gurobi.Optimizer(GUROBI_ENV, OutputFlag=0, QCPDual=0)
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
        @test_throws Gurobi.GurobiError MOI.get(model, MOI.ConstraintDual(), c1)
        @test_throws Gurobi.GurobiError MOI.get(model, MOI.ConstraintDual(), c2)
        @test_throws Gurobi.GurobiError MOI.get(model, MOI.ConstraintDual(), c3)
    end
end

@testset "Add constraints" begin
    model = Gurobi.Optimizer(GUROBI_ENV)
    x = MOI.add_variables(model, 2)
    MOI.add_constraints(
        model,
        [MOI.ScalarAffineFunction([MOI.ScalarAffineTerm(1.0, x[i])], 0.0) for i in 1:2],
        MOI.EqualTo.([0.0, 0.0])
    )
    @test MOI.get(model, MOI.NumberOfConstraints{
        MOI.ScalarAffineFunction{Float64}, MOI.EqualTo{Float64}
    }()) == 2
end

@testset "Extra name tests" begin
    model = Gurobi.Optimizer(GUROBI_ENV)
    @testset "Variables" begin
        MOI.empty!(model)
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

    @testset "Variable bounds" begin
        MOI.empty!(model)
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

    @testset "Affine constraints" begin
        MOI.empty!(model)
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
end

@testset "ConstraintAttribute" begin
    model = Gurobi.Optimizer(GUROBI_ENV)
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
    @test_throws MOI.UnsupportedAttribute(attr) MOI.set(model, attr, c2, 1)
    @test_throws MOI.UnsupportedAttribute(attr) MOI.get(model, attr, c2)
    # Setting an attribute to a value of the wrong type.
    @test_throws(
        ArgumentError("Attribute Lazy is Integer but Float64 provided."),
        MOI.set(model, Gurobi.ConstraintAttribute("Lazy"), c2, 1.0)
    )
end

@testset "VariableAttribute" begin
    model = Gurobi.Optimizer(GUROBI_ENV)
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
    @test_throws MOI.UnsupportedAttribute(attr) MOI.set(model, attr, x, 1)
    @test_throws MOI.UnsupportedAttribute(attr) MOI.get(model, attr, x)
    # Setting an attribute to a value of the wrong type.
    @test_throws(
        ArgumentError("Attribute BranchPriority is Integer but Float64 provided."),
        MOI.set(model, Gurobi.VariableAttribute("BranchPriority"), x, 1.0)
    )
end
