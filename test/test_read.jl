# Test MIP
#
#   maximize z
#
#   s.t.  3 p1 + 5 p2 + 4 p3 - z = 0
#         0.5 p1 + 2 p2 + 1 p3 <= 6
#
#         p1 is integer: 0 <= p1
#         p1 is integer: 0 <= p2
#         p1 is integer: 0 <= p3
#         z is binary

using Gurobi, Base.Test

@testset "GRBread" begin

    ## BUILD MODEL
    #-------------

    simple_model_env = Gurobi.Env()

    simple_model = Gurobi.Model(simple_model_env, "simple_mip", :maximize)

    add_ivar!(simple_model, 0., 0, Inf)  # p1
    add_ivar!(simple_model, 0., 0, Inf)  # p2
    add_ivar!(simple_model, 0., 0, Inf)  # p3
    add_cvar!(simple_model, 1., 0., Inf) # z
    update_model!(simple_model)

    add_constr!(simple_model, [3., 5., 4., -1.], '<', 0.)
    add_constr!(simple_model, [3., 5., 4., -1.], '>', 0.)

    add_constr!(simple_model, [0.5, 2., 1., 0.], '<', 6.)

    setparam!(simple_model, "Heuristics", 0.0)
    setparam!(simple_model, "Presolve", 0)

    update_model!(simple_model)

    optimize(simple_model)

    ## WRITE SOLUTION FILES
    #----------------------

    write_model(simple_model, "simple_out.sol")
    write_model(simple_model, "simple_out.mst")

    @test isfile("simple_out.sol")
    @test isfile("simple_out.mst")

    ## READ SOLUTION FILE
    #--------------------

    start_1 = Gurobi.get_dblattrarray(simple_model, "Start", 1, 4)
    @test start_1 == [1.0e101, 1.0e101, 1.0e101, 1.0e101]

    Gurobi.read(simple_model, "simple_out.sol")
    update_model!(simple_model)

    start_2 = Gurobi.get_dblattrarray(simple_model, "Start", 1, 4)
    @test start_2 == [12.0, 0.0, 0.0, 36.0]

    ## VERIFY MIP START
    #------------------

    #The parameters Presolve and Heuristics must be set to 0
    optimize(simple_model)

    ## DELETE FILES
    #--------------

    rm("simple_out.sol")
    rm("simple_out.mst")

end
