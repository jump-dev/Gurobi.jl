# QCQP example
#    maximize x + y
#
#    s.t.  x, y >= 0
#          x^2 + y^2 <= 1
#
#    solution: (0.71, 0.71) objv = 1.414

using Gurobi, Test

@testset "QCQP1" begin

    env = Gurobi.Env()
    setparam!(env, "OutputFlag", 0)

    model = Gurobi.Model(env, "qcqp_01", :maximize)

    add_cvars!(model, [1., 1.], 0., Inf)
    update_model!(model)

     # add_qpterms!(model, linearindices, linearcoeffs, qrowinds, qcolinds, qcoeffs, sense, rhs)
    add_qconstr!(model, [], [], [1, 2], [1, 2], [1, 1.], '<', 1.0)
    update_model!(model)

    optimize(model)

    @test isapprox(get_solution(model), [0.707107, 0.707107], atol=1e-5)
    @test isapprox(get_objval(model), 1.414214, atol=1e-5)
end
