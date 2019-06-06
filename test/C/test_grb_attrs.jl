# an example on mixed integer programming
#
#   maximize x + 2 y + 5 z
#
#   s.t.  x + y + z <= 10
#         x + 2 y + z <= 15
#
#         x is continuous: 0 <= x <= 5
#         y is integer: 0 <= y <= 10
#         z is binary
#

using Gurobi, Test

@testset "GRB Attributes" begin
    env = Gurobi.Env()

    model = Gurobi.Model(env, "mip_01", :maximize)

    add_cvar!(model, 1., 0., 5.)  # x
    add_ivar!(model, 2., 0, 10)   # y
    add_bvar!(model, 5.)          # z
    update_model!(model)

    add_constr!(model, ones(3), '<', 10.)
    add_constr!(model, [1., 2., 1.], '<', 15.)

    update_model!(model)
    @test Gurobi.get_intattr(model, "NumVars") == 3
    @test Gurobi.get_dblattr(model, "ObjCon") == 0.0
    @test Gurobi.get_strattr(model, "ModelName") == "mip_01"


    @test Gurobi.get_dblattrelement(model, "UB", 1) == 5.0
    @test Gurobi.get_charattrelement(model, "Sense", 1) == '<'

    Gurobi.set_dblattrelement!(model, "UB", 1, 4.0)
    Gurobi.set_charattrelement!(model, "Sense", 1, '>')
    update_model!(model)
    @test Gurobi.get_dblattrelement(model, "UB", 1) == 4.0
    @test Gurobi.get_charattrelement(model, "Sense", 1) == '>'
end
