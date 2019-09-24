using Gurobi, Test, Random

const MOI = Gurobi.MOI
const GUROBI_ENV = Gurobi.Env()

@testset "Lazy cut" begin
    model = Gurobi.Optimizer(
        GUROBI_ENV,
        OutputFlag = 0,
        Cuts = 0,
        Presolve = 0,
        Heuristics = 0,
        LazyConstraints = 1
    )
    MOI.Utilities.loadfromstring!(model, """
        variables: x, y
        maxobjective: y
        c1: x in Integer()
        c2: y in Integer()
        c3: x in Interval(0.0, 2.5)
        c4: y in Interval(0.0, 2.5)
    """)
    x = MOI.get(model, MOI.VariableIndex, "x")
    y = MOI.get(model, MOI.VariableIndex, "y")
    lazy_called = false
    MOI.set(model, MOI.LazyConstraintCallback(), cb_data -> begin
        lazy_called = true
        x_val = MOI.get(model, MOI.CallbackVariablePrimal(cb_data), x)
        y_val = MOI.get(model, MOI.CallbackVariablePrimal(cb_data), y)
        if y_val - x_val > 1 + 1e-6
            MOI.submit(
                model,
                MOI.LazyConstraint(cb_data),
                MOI.ScalarAffineFunction{Float64}(
                    MOI.ScalarAffineTerm.([-1.0, 1.0], [x, y]),
                    0.0
                ),
                MOI.LessThan{Float64}(1.0)
            )
        elseif y_val + x_val > 3 + 1e-6
            MOI.submit(
                model,
                MOI.LazyConstraint(cb_data),
                MOI.ScalarAffineFunction{Float64}(
                    MOI.ScalarAffineTerm.([1.0, 1.0], [x, y]),
                    0.0
                ), MOI.LessThan{Float64}(3.0)
            )
        end
    end)
    MOI.optimize!(model)
    @test lazy_called
    @test MOI.get(model, MOI.VariablePrimal(), x) == 1
    @test MOI.get(model, MOI.VariablePrimal(), y) == 2
end

@testset "User cut" begin
    model = Gurobi.Optimizer(
        GUROBI_ENV,
        OutputFlag = 0,
        Cuts = 0,
        Presolve = 0,
        PreCrush = 1,
        Heuristics = 0
    )
    N = 30
    x = MOI.add_variables(model, N)
    MOI.add_constraints(model, MOI.SingleVariable.(x), MOI.ZeroOne())
    MOI.set.(model, MOI.VariablePrimalStart(), x, 0.0)
    Random.seed!(1)
    item_weights, item_values = rand(N), rand(N)
    MOI.add_constraint(
        model,
        MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(item_weights, x), 0.0),
        MOI.LessThan(10.0)
    )
    MOI.set(
        model,
        MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(),
        MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(item_values, x), 0.0)
    )
    MOI.set(model, MOI.ObjectiveSense(), MOI.MAX_SENSE)

    user_cut_submitted = false
    MOI.set(model, MOI.UserCutCallback(), cb_data -> begin
        terms = MOI.ScalarAffineTerm{Float64}[]
        accumulated = 0.0
        for (i, xi) in enumerate(x)
            if MOI.get(model, MOI.CallbackVariablePrimal(cb_data), xi) > 0.0
                push!(terms, MOI.ScalarAffineTerm(1.0, xi))
                accumulated += item_weights[i]
            end
        end
        if accumulated > 10.0
            MOI.submit(
                model,
                MOI.UserCut(cb_data),
                MOI.ScalarAffineFunction{Float64}(terms, 0.0),
                MOI.LessThan{Float64}(length(terms) - 1)
            )
            user_cut_submitted = true
        end
    end)
    MOI.optimize!(model)
    @test user_cut_submitted
end

@testset "Heuristic Solution" begin
    model = Gurobi.Optimizer(
        GUROBI_ENV,
        OutputFlag = 0,
        Cuts = 0,
        Presolve = 0,
        Heuristics = 0
    )
    N = 30
    x = MOI.add_variables(model, N)
    MOI.add_constraints(model, MOI.SingleVariable.(x), MOI.ZeroOne())
    MOI.set.(model, MOI.VariablePrimalStart(), x, 0.0)
    Random.seed!(1)
    item_weights, item_values = rand(N), rand(N)
    MOI.add_constraint(
        model,
        MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(item_weights, x), 0.0),
        MOI.LessThan(10.0)
    )
    MOI.set(
        model,
        MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(),
        MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(item_values, x), 0.0)
    )
    MOI.set(model, MOI.ObjectiveSense(), MOI.MAX_SENSE)
    solution_accepted = false
    solution_rejected = false
    MOI.set(model, MOI.HeuristicCallback(), cb_data -> begin
        x_vals = MOI.get.(model, MOI.CallbackVariablePrimal(cb_data), x)
        if MOI.submit(
            model,
            MOI.HeuristicSolution(cb_data),
            x,
            floor.(x_vals)
        ) == MOI.HEURISTIC_SOLUTION_ACCEPTED
            solution_accepted = true
        end
        if MOI.submit(
            model,
            MOI.HeuristicSolution(cb_data),
            x,
            ceil.(x_vals)
        ) == MOI.HEURISTIC_SOLUTION_REJECTED
            solution_rejected = true
        end
    end)
    MOI.optimize!(model)
    @test solution_accepted
    @test solution_rejected
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
            return nothing
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
                Gurobi.cbget_mipsol_sol(m, cb_data, cb_where)
                x_val = MOI.get(m, MOI.CallbackVariablePrimal(cb_data), x)
                y_val = MOI.get(m, MOI.CallbackVariablePrimal(cb_data), y)
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
                    MOI.submit(m, MOI.LazyConstraint(cb_data),
                        MOI.ScalarAffineFunction{Float64}(
                            MOI.ScalarAffineTerm.([-1.0, 1.0], [x, y]),
                            0.0
                        ),
                        MOI.LessThan{Float64}(1.0)
                    )
                elseif y_val + x_val > 3 + TOL
                    MOI.submit(m, MOI.LazyConstraint(cb_data),
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
