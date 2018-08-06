using Gurobi, Base.Test, MathOptInterface, MathOptInterface.Test

const MOI  = MathOptInterface
const MOIT = MathOptInterface.Test

@testset "Unit Tests" begin
    config = MOIT.TestConfig()
    solver = GurobiOptimizer(OutputFlag=0)

    # TODO(@odow): see MathOptInterface Issue #404
    # The basic constraint tests incorrectly add multiple constraints
    # that are illegal, e.g., two SingleVariable-in-ZeroOne constraints
    # for the same variable.
    # MOIT.basic_constraint_tests(solver, config)

    MOIT.unittest(solver, config, [
        "solve_affine_interval",  # Interval constraints not wrapped
        "solve_qcp_edge_cases"    # tested below
    ])

    @testset "solve_qcp_edge_cases" begin
        MOIT.solve_qcp_edge_cases(solver,
            MOIT.TestConfig(atol=1e-4)
        )
    end

    MOIT.modificationtest(solver, config, [
        "solve_func_scalaraffine_lessthan"
    ])
end

@testset "Linear tests" begin
    linconfig = MOIT.TestConfig()
    @testset "Default Solver"  begin
        solver = GurobiOptimizer(OutputFlag=0)
        MOIT.contlineartest(solver, linconfig, ["linear10","linear12","linear8a","linear8b","linear8c"])
    end
    @testset "InfUnbdInfo=1" begin
        solver_nopresolve = GurobiOptimizer(OutputFlag=0, InfUnbdInfo=1)
        MOIT.contlineartest(solver_nopresolve, linconfig, ["linear10","linear12","linear8a"])
    end
    @testset "No certificate" begin
        solver = GurobiOptimizer(OutputFlag=0)
        linconfig_nocertificate = MOIT.TestConfig(infeas_certificates=false)
        MOIT.linear12test(solver, linconfig_nocertificate)
        MOIT.linear8atest(solver, linconfig_nocertificate)
    end
    # 10 is ranged
end

@testset "Quadratic tests" begin
    quadconfig = MOIT.TestConfig(atol=1e-4, rtol=1e-4, duals=false, query=false)
    solver = GurobiOptimizer(OutputFlag=0)
    MOIT.contquadratictest(solver, quadconfig)
end

@testset "Linear Conic tests" begin
    linconfig = MOIT.TestConfig()
    solver = GurobiOptimizer(OutputFlag=0)
    MOIT.lintest(solver, linconfig, ["lin3","lin4"])

    solver_nopresolve = GurobiOptimizer(OutputFlag=0, InfUnbdInfo=1)
    MOIT.lintest(solver_nopresolve, linconfig)
end

@testset "Integer Linear tests" begin
    intconfig = MOIT.TestConfig()
    solver = GurobiOptimizer(OutputFlag=0)
    MOIT.intlineartest(solver, intconfig, ["int3"])

    # 3 is ranged
end
@testset "ModelLike tests" begin
    solver = GurobiOptimizer()
    MOIT.nametest(solver)
    @testset "validtest" begin
        MOIT.validtest(solver)
    end
    @testset "emptytest" begin
        MOIT.emptytest(solver)
    end
    @testset "orderedindicestest" begin
        # TODO(@odow): see MathOptInterface Issue #404
        # The basic constraint tests incorrectly add multiple constraints
        # that are illegal, e.g., two SingleVariable-in-ZeroOne constraints
        # for the same variable.
        # MOIT.orderedindicestest(solver)
    end
    @testset "canaddconstrainttest" begin
        MOIT.canaddconstrainttest(solver, Float64, Complex{Float64})
    end
    @testset "copytest" begin
        solver2 = GurobiOptimizer()
        MOIT.copytest(solver,solver2)
    end
end

@testset "Gurobi Callback" begin
    @testset "Generic callback" begin
        m = GurobiOptimizer(OutputFlag=0)
        x = MOI.addvariable!(m)
        MOI.addconstraint!(m, MOI.SingleVariable(x), MOI.GreaterThan(1.0))
        MOI.set!(m, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(),
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

        MOI.set!(m, Gurobi.CallbackFunction(), callback_function)
        MOI.optimize!(m)

        @test length(cb_calls) > 0
        @test Gurobi.CB_MESSAGE in cb_calls
        @test Gurobi.CB_PRESOLVE in cb_calls
        @test !(Gurobi.CB_MIPSOL in cb_calls)
    end

    @testset "Lazy cut" begin
        m = GurobiOptimizer(OutputFlag=0, Cuts=0, Presolve=0, Heuristics=0, LazyConstraints=1)
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
                Gurobi.loadcbsolution!(m, cb_data, cb_where)
                x_val = MOI.get(m, MOI.VariablePrimal(), x)
                y_val = MOI.get(m, MOI.VariablePrimal(), y)
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

        MOI.set!(m, Gurobi.CallbackFunction(), callback_function)
        MOI.optimize!(m)

        @test MOI.get(m, MOI.VariablePrimal(), x) == 1
        @test MOI.get(m, MOI.VariablePrimal(), y) == 2

        @test length(cb_calls) > 0
        @test Gurobi.CB_MESSAGE in cb_calls
        @test Gurobi.CB_PRESOLVE in cb_calls
        @test Gurobi.CB_MIPSOL in cb_calls
    end
end

@testset "User limit handling (issue #140)" begin
    # Verify that we return the correct status codes when a mixed-integer
    # problem has been solved to a *feasible* but not necessarily optimal
    # solution. To do that, we will set up an intentionally dumbed-down
    # Gurobi optimizer (with all heuristics and pre-solve turned off) and
    # ask it to solve a classic knapsack problem. Setting SolutionLimit=1
    # forces the solver to return after its first feasible MIP solution,
    # which tests the right part of the code without relying on potentially
    # flaky or system-dependent time limits.
    m = GurobiOptimizer(OutputFlag=0,
                        SolutionLimit=1,
                        Heuristics=0.0,
                        Presolve=0)
    N = 100
    x = MOI.addvariables!(m, N)
    for xi in x
        MOI.addconstraint!(m, MOI.SingleVariable(xi), MOI.ZeroOne())
        MOI.set!(m, MOI.VariablePrimalStart(), xi, 0.0)
    end
    # Given a collection of items with individual weights and values,
    # maximize the total value carried subject to the constraint that
    # the total weight carried is less than 10.
    srand(1)
    item_weights = rand(N)
    item_values = rand(N)
    MOI.addconstraint!(m,
        MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(item_weights, x), 0.0),
        MOI.LessThan(10.0))
    MOI.set!(m, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(),
        MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(.-item_values, x), 0.0))
    MOI.optimize!(m)

    @test MOI.get(m, MOI.TerminationStatus()) == MOI.SolutionLimit
    # We should have a primal feasible solution:
    @test MOI.get(m, MOI.PrimalStatus()) == MOI.FeasiblePoint
    # But we have no dual status:
    @test MOI.get(m, MOI.DualStatus()) == MOI.UnknownResultStatus
end

@testset "Constant objective (issue #111)" begin
    m = GurobiOptimizer()
    x = MOI.addvariable!(m)
    MOI.set!(m, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(),
        MOI.ScalarAffineFunction(MOI.ScalarAffineTerm{Float64}[], 2.0))
    @test MOI.get(m, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}()).constant == 2.0
    @test Gurobi.get_dblattr(m.inner, "ObjCon") == 2.0

    MOI.modify!(m, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(), MOI.ScalarConstantChange(3.0))
    @test MOI.get(m, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}()).constant == 3.0
    @test Gurobi.get_dblattr(m.inner, "ObjCon") == 3.0
end
