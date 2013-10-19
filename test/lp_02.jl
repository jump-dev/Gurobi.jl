# LP programming using MATLAB-like construction
#
#   maximize 1000 x + 350 y
#
#       s.t. x >= 30, y >= 0
#            x - 1.5y >= 0  (i.e. -x + 1.5 y <= 0)
#            12 x + 8 y <= 1000
#            1000 x + 300 y <= 70000
#            
#   solution: (59.0909, 36.3636)
#   objv: 71818.1818
#

using MathProgBase
using Gurobi

env = Gurobi.Env()

model = gurobi_model(env; 
	name="lp_02", 
	sense=:maximize, 
	f = [1000., 350.],
	A = [-1. 1.5; 12. 8.; 1000. 300.], 
	b = [0., 1000., 70000.], 
	lb = [0., 30.])

println(model)

optimize(model)

println()
println("soln = $(get_solution(model))")
println("objv = $(get_objval(model))")
