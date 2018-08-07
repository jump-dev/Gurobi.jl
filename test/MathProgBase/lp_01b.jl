# a simple LP example
#
#   maximize x + y
#
#   s.t. 50 x + 24 y <= 2400
#        30 x + 33 y <= 2100
#        x >= 45, y >= 5
#
#   solution: x = 45, y = 6.25, objv = 51.25

using Gurobi, Compat.Test, Compat.GC

@testset "LP 01b" begin

    env = Gurobi.Env()
    setparam!(env, "OutputFlag", 0)
    setparams!(env, Method=2)  # using barrier method

    method = getparam(env, "Method")

    model = Gurobi.Model(env, "lp_01", :maximize)

    # add variables
    add_cvars!(model, [1., 1.], [45., 5.], Inf)
    update_model!(model)

    # add constraints
    add_constrs!(model, Cint[1, 3], Cint[1, 2, 1, 2],
        [50., 24., 30., 33.], '<', [2400., 2100.])
    update_model!(model)

    # perform optimization
    optimize(model)

    # show results
    info = get_optiminfo(model)

    sol = get_solution(model)
    @test sol == [45, 6.25]

    objv = get_objval(model)
    @test objv == 51.25

    GC.gc()  # test finalizers
end
