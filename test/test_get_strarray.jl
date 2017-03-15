using Gurobi
using Base.Test

# enviroment and problem
env = Gurobi.Env()
m = Gurobi.Model(env, "test")

# read mps problem in Gurobi
read_model(m, "test_get_strarray.mps")

# test functions
@test Gurobi.get_strattrarray(m, "VarName", 1, Gurobi.num_vars(m))[1] == "XONE"
@test Gurobi.get_strattrarray(m, "ConstrName", 1, Gurobi.num_constrs(m))[1] =="LIM1"

# get variable and constraint names
varnames = Gurobi.get_strattrarray(m, "VarName", 1, Gurobi.num_vars(m))
connames = Gurobi.get_strattrarray(m, "ConstrName", 1, Gurobi.num_constrs(m))

# print variable and constraint names
println(varnames)		# should output ["XONE","YTWO","ZTHREE"]
println(connames)		# should output ["LIM1","LIM2","MYEQN"]
