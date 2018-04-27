# a simple LP example
#
#   maximize x + y
#
#   s.t. 50 x + 24 y <= 2400
#        30 x + 33 y <= 2100
#        x >= 45, y >= 5
#
#   solution: x = 45, y = 6.25, objv = 51.25

using Gurobi, Base.Test

@testset "LP 01a" begin
    env = Gurobi.Env()

    method = getparam(env, "Method")
    println("method = $method")

    model = Gurobi.Model(env, "lp_01", :maximize)

    # add variables
    add_cvar!(model, 1.0, 45., Inf)  # x
    add_cvar!(model, 1.0,  5., Inf)  # y
    update_model!(model)

    # add constraints
    add_constr!(model, [50., 24.], '<', 2400.)
    add_constr!(model, [30., 33.], '<', 2100.)
    update_model!(model)

    println(model)

    # perform optimization
    optimize(model)

    # show results
    info = get_optiminfo(model)
    println()
    println(info)

    @test get_solution(model) == [45, 6.25]
    @test get_objval(model) == 51.25

    gc()  # test finalizers
end
