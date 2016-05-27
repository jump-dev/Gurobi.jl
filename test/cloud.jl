using Base.Test
using Gurobi
using MathProgBase.SolverInterface

env = Gurobi.Env("csdemo.gurobi.com")
Gurobi.free_env(env)

@test_throws ErrorException Gurobi.Env("csdemo.gurobi.com", "wrongpassword")

@test_throws ErrorException Gurobi.Env("fakehostdoesnotexist", "wrongpassword")

solver = GurobiSolver(cloudhost="csdemo.gurobi.com")
LinearQuadraticModel(solver)
