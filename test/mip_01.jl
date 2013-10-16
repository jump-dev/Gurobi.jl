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

using Gurobi

env = Gurobi.Env()

model = Gurobi.Model(env, "mip_01", :maximize)

add_cvar!(model, 1., 0., 5.)  # x
add_ivar!(model, 2., 0, 10)   # y
add_bvar!(model, 5.)          # z
update_model!(model)

add_constr!(model, ones(3), '<', 10.)
add_constr!(model, [1., 2., 1.], '<', 15.)

println(model)

optimize(model)

println("sol = $(get_solution(model))")
println("objv = $(get_objval(model))")
