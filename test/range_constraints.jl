using Gurobi, MathProgBase

# Test Range constraints pass
m = MathProgBase.LinearQuadraticModel(GurobiSolver())
MathProgBase.loadproblem!(m, [1 1], [0, 0],[1,1], [1,1], [0], [1], :Max)
MathProgBase.optimize!(m)
