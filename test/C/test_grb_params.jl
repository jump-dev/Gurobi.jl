using Test, Gurobi

@testset "GRB Parameters" begin
    env = Gurobi.Env()
    model = Gurobi.Model(env, "")

    @test_throws Gurobi.GurobiError Gurobi.get_str_param(model, "Foobar")
    @test_throws Gurobi.GurobiError Gurobi.set_str_param!(model, "Foobar", "")
    @test_throws Gurobi.GurobiError Gurobi.get_int_param(model, "Foobar")
    @test_throws Gurobi.GurobiError Gurobi.set_int_param!(model, "Foobar", 0)
    @test_throws Gurobi.GurobiError Gurobi.get_dbl_param(model, "Foobar")
    @test_throws Gurobi.GurobiError Gurobi.set_dbl_param!(model, "Foobar", 0.0)
end
