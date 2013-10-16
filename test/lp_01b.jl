# a simple LP example
#
#   maximize x + y
#
#   s.t. 50 x + 24 y <= 2400
#        30 x + 33 y <= 2100
#        x >= 45, y >= 5
#
#   solution: x = 45, y = 6.25, objv = 51.25

using Gurobi

env = Gurobi.Env()
set_int_param!(env, "Method", 2)  # using barrier method

model = Gurobi.Model(env, "lp_01", :maximize)

# add variables
add_cvars!(model, [1., 1.], [45., 5.], nothing)
update_model!(model)

# add constraints
add_constrs!(model, Cint[1, 3], Cint[1, 2, 1, 2], 
    [50., 24., 30., 33.], '<', [2400., 2100.])
update_model!(model)

println(model)

# perform optimization
optimize(model)

# show results
info = get_optim_info(model)
println()
println(info)

sol = get_solution(model)
println("soln = $(sol)")

objv = get_objval(model)
println("objv = $(objv)")

gc()  # test finalizers
