module TestMultiobjective

using Gurobi
using Test

function runtests()
    for name in names(@__MODULE__; all = true)
        if startswith("$(name)", "test_")
            @testset "$(name)" begin
                getfield(@__MODULE__, name)()
            end
        end
    end
    return
end

const MOI = Gurobi.MOI

const GRB_ENV = isdefined(Main, :GRB_ENV) ? Main.GRB_ENV : Gurobi.Env()

function test_multiobjective()
    model = Gurobi.Optimizer(GRB_ENV)
    MOI.set(model, MOI.Silent(), true)
    MOI.Utilities.loadfromstring!(
        model,
        """
variables: x, y
minobjective: 2x + y
c1: x + y >= 1.0
c2: 0.5 * x + 1.0 * y >= 0.75
c3: x >= 0.0
c4: y >= 0.25
""",
    )
    x = MOI.get(model, MOI.VariableIndex, "x")
    y = MOI.get(model, MOI.VariableIndex, "y")

    f = MOI.ScalarAffineFunction(
        [MOI.ScalarAffineTerm(1.0, x), MOI.ScalarAffineTerm(3.0, y)],
        0.0,
    )

    MOI.set(model, Gurobi.MultiObjectiveFunction(2), f)

    @test MOI.get(model, Gurobi.MultiObjectiveWeight(1)) == 1.0
    @test MOI.get(model, Gurobi.MultiObjectiveWeight(2)) == 1.0
    @test MOI.get(model, Gurobi.MultiObjectivePriority(1)) == 0
    @test MOI.get(model, Gurobi.MultiObjectivePriority(2)) == 0

    MOI.optimize!(model)
    @test MOI.get(model, MOI.ObjectiveValue()) ≈ 1.5
    @test MOI.get(model, MOI.VariablePrimal(), x) ≈ 0.5
    @test MOI.get(model, MOI.VariablePrimal(), y) ≈ 0.5

    BFS = [
        (x = 1.0, y = 0.25, f1 = 2.25, f2 = 1.75),
        (x = 0.5, y = 0.5, f1 = 1.5, f2 = 2.0),
        (x = 0.0, y = 1.0, f1 = 1.0, f2 = 3.0),
    ]
    for (i, λ) in enumerate([0.2, 0.5, 0.8])
        MOI.set(model, Gurobi.MultiObjectiveWeight(1), λ)
        MOI.set(model, Gurobi.MultiObjectiveWeight(2), 1 - λ)
        MOI.optimize!(model)
        @test MOI.get(model, MOI.VariablePrimal(), x) ≈ BFS[i].x
        @test MOI.get(model, MOI.VariablePrimal(), y) ≈ BFS[i].y
        @test MOI.get(model, Gurobi.MultiObjectiveValue(1)) ≈ BFS[i].f1
        @test MOI.get(model, Gurobi.MultiObjectiveValue(2)) ≈ BFS[i].f2
    end

    MOI.set(model, Gurobi.MultiObjectiveWeight(1), 1.0)
    MOI.set(model, Gurobi.MultiObjectiveWeight(2), 1.0)
    MOI.set(model, Gurobi.MultiObjectivePriority(1), 1)
    MOI.set(model, Gurobi.MultiObjectivePriority(2), 2)
    @test MOI.get(model, Gurobi.MultiObjectivePriority(1)) == 1
    @test MOI.get(model, Gurobi.MultiObjectivePriority(2)) == 2

    MOI.optimize!(model)

    @test MOI.get(model, MOI.VariablePrimal(), x) ≈ BFS[1].x
    @test MOI.get(model, MOI.VariablePrimal(), y) ≈ BFS[1].y
    @test MOI.get(model, Gurobi.MultiObjectiveValue(1)) ≈ BFS[1].f1
    @test MOI.get(model, Gurobi.MultiObjectiveValue(2)) ≈ BFS[1].f2

    MOI.set(model, Gurobi.MultiObjectivePriority(1), 2)
    MOI.set(model, Gurobi.MultiObjectivePriority(2), 1)
    @test MOI.get(model, Gurobi.MultiObjectivePriority(1)) == 2
    @test MOI.get(model, Gurobi.MultiObjectivePriority(2)) == 1

    MOI.optimize!(model)

    @test MOI.get(model, MOI.VariablePrimal(), x) ≈ BFS[3].x
    @test MOI.get(model, MOI.VariablePrimal(), y) ≈ BFS[3].y
    @test MOI.get(model, Gurobi.MultiObjectiveValue(1)) ≈ BFS[3].f1
    @test MOI.get(model, Gurobi.MultiObjectiveValue(2)) ≈ BFS[3].f2
end

end

TestMultiobjective.runtests()
