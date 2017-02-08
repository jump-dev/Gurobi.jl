# hierarchical multi-objetive LP
# 

using Gurobi

env = Gurobi.Env()

A = [1. 1.]
lb = [0.; 0.]
b = [1.]
f = [1.; 1.]
model = gurobi_model(env, sense = :maximize, A = A, b = b, lb = lb, f = [0., 0.])

set_multiobj_n!(model, 2)
set_multiobj!(model, 0, [1., 1.], 10, 1.0)

set_multiobj!(model, 1, [1., 0.], 1, 1.0)
update_model!(model)

optimize(model)
Test.@test get_solution(model) == [1., 0.]

set_multiobj!(model, 1, [0., 1.], 1, 1.0)
update_model!(model)

optimize(model)
Test.@test get_solution(model) == [0., 1.]
