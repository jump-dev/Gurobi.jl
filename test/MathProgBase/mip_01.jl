# an example on mixed integer programming
#
#   maximize x + 3 y + 5 z
#
#   s.t.  x + y + z <= 10
#         x + 2 y + z <= 15
#
#         x is continuous: 0 <= x <= 5
#         y is integer: 0 <= y <= 10
#         z is binary
#

using Gurobi, Compat.Test

@testset "MIP 01" begin
    env = Gurobi.Env()
    setparam!(env, "OutputFlag", 0)

    model = Gurobi.Model(env, "mip_01", :maximize)

    add_cvar!(model, 1., 0., 5.)  # x
    add_ivar!(model, 3., 0, 10)   # y
    add_bvar!(model, 5.)          # z
    update_model!(model)

    add_constr!(model, ones(3), '<', 10.)
    add_constr!(model, [1., 2., 1.], '<', 15.)

    optimize(model)
    @test get_solution(model) == [0.0, 7.0, 1.0]
    @test get_objval(model) == 26.0
end
