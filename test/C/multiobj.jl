# hierarchical multi-objetive LP
#

using Gurobi, Compat.Test

env = Gurobi.Env()
setparam!(env, "OutputFlag", 0)

A = [1. 1.]
lb = [0.; 0.]
b = [1.]
f = [1.; 1.]
model = gurobi_model(env, sense = :maximize, A = A, b = b, lb = lb, f = [0., 0.])

Gurobi.set_multiobj_n!(model, 2)
Gurobi.update_model!(model)
Gurobi.set_multiobj!(model, 1, [1., 1.], 10, 1.0)

Gurobi.set_multiobj!(model, 2, [1., 0.], 1, 1.0)
Gurobi.update_model!(model)

optimize(model)
@test get_solution(model) == [1., 0.]

Gurobi.set_multiobj!(model, 2, [0., 1.], 1, 1.0)
Gurobi.update_model!(model)

optimize(model)
@test get_solution(model) == [0., 1.]
