# QP example
#
#    minimize 2 x^2 + y^2 + xy + x + y
#
#       s.t.  x, y >= 0
#             x + y = 1
#
#    solution: (0.25, 0.75), objv = 1.875
#

using Gurobi, Compat.Test

@testset "QP1" begin

    env = Gurobi.Env()
    setparam!(env, "OutputFlag", 0)

    model = Gurobi.Model(env, "qp_01")

    add_cvars!(model, [1., 1.], 0., Inf)
    update_model!(model)

    add_qpterms!(model, [1, 1, 2], [1, 2, 2], [2., 1., 1.])
    add_constr!(model, [1., 1.], '=', 1.)
    update_model!(model)

    optimize(model)

    @test isapprox(get_solution(model), [0.25, 0.75], atol=1e-4)
    @test get_objval(model) ≈ 1.875
end
