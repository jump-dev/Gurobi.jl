# Copyright (c) 2015 Dahua Lin, Miles Lubin, Joey Huchette, Iain Dunning, and
#   contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

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

const GRB_ENV =
    isdefined(Main, :GRB_ENV) ? Main.GRB_ENV : Gurobi.Env(output_flag = 0)

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
    @test MOI.get(model, MOI.ObjectiveValue()) ≈ [1.5, 2.0]
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
    MOI.set(model, Gurobi.MultiObjectiveWeight(1), 1)
    MOI.set(model, Gurobi.MultiObjectiveWeight(2), 1)
    @test MOI.get(model, Gurobi.MultiObjectiveWeight(1)) == 1.0
    @test MOI.get(model, Gurobi.MultiObjectiveWeight(2)) == 1.0
    return
end

function test_example_biobjective_knapsack()
    p1 = [77.0, 94, 71, 63, 96, 82, 85, 75, 72, 91, 99, 63, 84, 87, 79, 94, 90]
    p2 = [65.0, 90, 90, 77, 95, 84, 70, 94, 66, 92, 74, 97, 60, 60, 65, 97, 93]
    w = [80.0, 87, 68, 72, 66, 77, 99, 85, 70, 93, 98, 72, 100, 89, 67, 86, 91]
    model = Gurobi.Optimizer()
    x = MOI.add_variables(model, length(w))
    MOI.add_constraint.(model, x, MOI.ZeroOne())
    MOI.add_constraint(model, w' * x, MOI.LessThan(900.0))
    obj_f = MOI.Utilities.operate(vcat, Float64, p1' * x, p2' * x)
    MOI.set(model, MOI.ObjectiveSense(), MOI.MAX_SENSE)
    MOI.set(model, MOI.ObjectiveFunction{typeof(obj_f)}(), obj_f)
    MOI.optimize!(model)
    results = Dict(
        [955.0, 906.0] => [2, 3, 5, 6, 9, 10, 11, 14, 15, 16, 17],
        [948.0, 939.0] => [1, 2, 3, 5, 6, 8, 10, 11, 15, 16, 17],
        [934.0, 971.0] => [2, 3, 5, 6, 8, 10, 11, 12, 15, 16, 17],
        [918.0, 983.0] => [2, 3, 4, 5, 6, 8, 10, 11, 12, 16, 17],
    )
    found_non_dominated_point = false
    for i in 1:MOI.get(model, MOI.ResultCount())
        X = findall(elt -> elt > 0.9, MOI.get.(model, MOI.VariablePrimal(i), x))
        Y = MOI.get(model, MOI.ObjectiveValue(i))
        if haskey(results, Y)
            @test results[Y] == X
            found_non_dominated_point = true
        end
    end
    @test found_non_dominated_point
    return
end

end

TestMultiobjective.runtests()
