#   maximize 4 x + 3 y
#
#   s.t. 2 x + 1 y <= 4
#        1 x + 2 y <= 4
#        x >= 0, y >= 0, z >= 0
#
#   solution: x = 1, y = 1.5, objv = 8.5

using Gurobi, Base.Test

@testset "GRB basis" begin
    env = Gurobi.Env()

    model = Gurobi.Model(env, "test basis", :maximize)

    add_cvar!(model, 4., 0., Inf)  # x
    add_cvar!(model, 3., 0., Inf)  # y
    update_model!(model)

    add_constr!(model, [2., 1., 0.], '<', 4.)
    add_constr!(model, [1., 2., 0.], '<', 1.)

    update_model!(model)

    println(model)

    optimize(model)
    @test get_solution(model) ≈ [1, 0.]
    @test get_objval(model) ≈ 4.

    @test Gurobi.get_basisidx(model) == [1, 3]
    @test get_tableaurow(model, 1) ≈ [1., 2., 0., 1.]
    @test get_tableaurow(model, 2) ≈ [0., -3., 1., -2.]
end
