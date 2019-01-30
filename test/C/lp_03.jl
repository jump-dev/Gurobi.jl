# Test get/set objective coefficients in LP

using Gurobi, Compat.Test, Compat.GC

@testset "LP 03" begin
	env = Gurobi.Env()
	setparam!(env, "OutputFlag", 0)

	# original model
	#
	#   maximize  2x + 2y
	#
	#	s.t. 0.2 <= x, y <= 1
	#

	model = gurobi_model(env;
		name="lp_03",
		sense=:maximize,
		f=[2.0, 2.0],
		lb=[0.2, 0.2],
		ub=[1.0, 1.0])

	lb_ = lowerbounds(model)
	ub_ = upperbounds(model)
	c_ = objcoeffs(model)

	@test lb_ == [0.2, 0.2]
	@test ub_ == [1.0, 1.0]
	@test c_ == [2.0, 2.0]

	optimize(model)

	@test get_solution(model) == [1.0, 1.0]
	@test get_objval(model) == 4.0

	# change objective (warm start)
	#
	#	maximize x - y
	#
	#	s.t. 0.2 <= x, y <= 1
	#

	set_objcoeffs!(model, [1, -1])
	update_model!(model)

	c_ = objcoeffs(model)
	@test c_ == [1.0, -1.0]

	optimize(model)

	@test get_solution(model) == [1.0, 0.2]
	@test get_objval(model) == 0.8

	GC.gc()
end
