# a simple LP example

using Gurobi

env = Gurobi.Env()

model = gurobi_model(env, "lp_01")
add_cvars!(model, [1., 2., 3.], [1., 2., 3.], nothing)
update_model!(model)

println(model)

optimize(model)

info = get_optim_info(model)
println(info)

sol = get_solution(model)
println("soln = $(sol)")

objv = get_objval(model)
println("objv = $(objv)")

gc()  # test finalizers
