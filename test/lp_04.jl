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

using Gurobi

env = Gurobi.Env()
setparams!(env, Method=2)  # using barrier method

method = getparam(env, "Method")
println("method = $method")

model = Gurobi.Model(env, "lp_01", :maximize)

# add variables
add_cvars!(model, [1., 1.], [0., 0.], Inf)
update_model!(model)

# add constraints
add_constrs!(model, Cint[1, 3], Cint[1, 2, 1, 2], 
    [2., 1., 1., 2.], '<', [4., 4.])
update_model!(model)

println(model)

# perform optimization
optimize(model)

# show results
info = get_optiminfo(model)
println()
println(info)

sol = get_solution(model)
println("soln = $(sol)")

objv = get_objval(model)
println("objv = $(objv)")

# PART 2: 
# copy and solve

model2 = copy(model)

println(model2)

# perform optimization
optimize(model2)

# show results
info = get_optiminfo(model2)
println()
println(info)

sol = get_solution(model2)
println("soln = $(sol)")

objv = get_objval(model2)
println("objv = $(objv)")

# PART 3: 
# change coeff and solve

#   maximize x + y
#
#   s.t. 2 x + 2 y <= 4
#        1 x + 2 y <= 4
#        x >= 0, y >= 0
#
#   solution: x = 0, y = 2, objv = 2

chg_coeffs!(model, [1], [2],  [2.])
update_model!(model)

println(model)

# perform optimization
optimize(model)

# show results
info = get_optiminfo(model)
println()
println(info)

sol = get_solution(model)
println("soln = $(sol)")

objv = get_objval(model)
println("objv = $(objv)")

# PART 4: 
# change coeff and solve

#   maximize x + y
#
#   s.t. 1 x + 2 y <= 4
#        x >= 0, y >= 0
#
#   solution: x = 4, y = 0, objv = 4

del_constrs!(model, [1])
update_model!(model)

println(model)

# perform optimization
optimize(model)

# show results
info = get_optiminfo(model)
println()
println(info)

sol = get_solution(model)
println("soln = $(sol)")

objv = get_objval(model)
println("objv = $(objv)")

# PART 5: 
# change coeff and solve

#   maximize y
#
#   s.t.  2 y <= 4
#           y >= 0
#
#   solution: y = 2, objv = 4

del_vars!(model, [1])
update_model!(model)

println(model)

# perform optimization
optimize(model)

# show results
info = get_optiminfo(model)
println()
println(info)

sol = get_solution(model)
println("soln = $(sol)")

objv = get_objval(model)
println("objv = $(objv)")



gc()  # test finalizers