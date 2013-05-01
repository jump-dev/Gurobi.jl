# QCQP example 
#    maximize x + y
#
#    s.t.  x, y >= 0
#          x^2 + y^2 <= 1
#
#    solution: (0.71, 0.71) objv = 1.414

using Gurobi

env = Gurobi.Env()

model = gurobi_model(env, "qcqp_01", :maximize)

add_cvars!(model, [1., 1.], 0., Inf)
update_model!(model)

 # add_qpterms!(model, linearindices, linearcoeffs, qrowinds, qcolinds, qcoeffs, sense, rhs)
add_qconstr!(model, [], [], [1, 2], [1, 2], [1, 1.], '<', 1.0)
update_model!(model)

println(model)

optimize(model)

println("sol = $(get_solution(model))")
println("obj = $(get_objval(model))")

