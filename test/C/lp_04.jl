# Start from simple LP
# Solve it
# Copy and solve again
# Chg coeff, solve, change back solve
# del constr and solve
# del var and solve

# a simple LP example
#
#   maximize x + y
#
#   s.t. 2 x + 1 y <= 4
#        1 x + 2 y <= 4
#        x >= 0, y >= 0
#
#   solution: x = 1.3333333, y = 1.3333333, objv = 2.66666666

using Gurobi, Compat.Test, Compat.GC

@testset "LP 04" begin
    env = Gurobi.Env()
    setparam!(env, "OutputFlag", 0)
    setparams!(env, Method=2)  # using barrier method

    method = getparam(env, "Method")
    @test method == 2

    model = Gurobi.Model(env, "lp_01", :maximize)

    # add variables
    add_cvars!(model, [1., 1.], [0., 0.], Inf)
    update_model!(model)

    # add constraints
    add_constrs!(model, Cint[1, 3], Cint[1, 2, 1, 2],
        [2., 1., 1., 2.], '<', [4., 4.])
    update_model!(model)

    # perform optimization
    optimize(model)

    @test isapprox(get_solution(model), [1.3333, 1.3333], atol=1e-4)
    @test isapprox(get_objval(model), 2.6666, atol=1e-4)

    # PART 2:
    # copy and solve

    model2 = copy(model)
    optimize(model2)
    @test isapprox(get_solution(model2), [1.3333, 1.3333], atol=1e-4)
    @test isapprox(get_objval(model2), 2.6666, atol=1e-4)

    # PART 3:
    # change coeff and solve

    #   maximize x + y
    #
    #   s.t. 2 x + 2 y <= 4
    #        1 x + 2 y <= 4
    #        x >= 0, y >= 0
    #
    #   solution: x = 0, y = 2, objv = 2
    #        or : any λ(0, 2) + (1-λ)(2, 0)

    chg_coeffs!(model, [1], [2],  [2.])
    update_model!(model)
    optimize(model)
    sol = get_solution(model)
    @test !any(sol .< 0)
    @test !any(sol .> 2)
    @test sum(sol) == 2.0
    @test isapprox(get_objval(model), 2.0, atol=1e-4)

    # PART 4:
    # del constr and solve

    #   maximize x + y
    #
    #   s.t. 1 x + 2 y <= 4
    #        x >= 0, y >= 0
    #
    #   solution: x = 4, y = 0, objv = 4

    del_constrs!(model, [1])
    update_model!(model)
    optimize(model)
    @test isapprox(get_solution(model), [4.0, 0.0], atol=1e-4)
    @test isapprox(get_objval(model), 4.0, atol=1e-4)

    # PART 5:
    # del var and solve

    #   maximize y
    #
    #   s.t.  2 y <= 4
    #           y >= 0
    #
    #   solution: y = 2, objv = 2

    del_vars!(model, [1])
    update_model!(model)
    optimize(model)
    @test isapprox(get_solution(model), [2.0], atol=1e-4)
    @test isapprox(get_objval(model), 2.0, atol=1e-4)

    GC.gc()  # test finalizers
end
