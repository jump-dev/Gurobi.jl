using Gurobi, Test

@testset "GRB String Attributes" begin
    # enviroment and problem
    env = Gurobi.Env()
    m = Gurobi.Model(env, "test")

    # read mps problem in Gurobi
    read_model(m, joinpath(@__DIR__, "test_get_strarray.mps"))

    # test functions
    @test Gurobi.get_strattrarray(m, "VarName", 1, Gurobi.num_vars(m)) == ["XONE","YTWO","ZTHREE"]
    @test Gurobi.get_strattrarray(m, "ConstrName", 1, Gurobi.num_constrs(m)) == ["LIM1","LIM2","MYEQN"]
end
