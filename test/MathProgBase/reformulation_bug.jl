using Gurobi, MathProgBase

# This problem should error out on Gurobi Versions 6.5.2 and earlier,
# but pass on later versions

m = MathProgBase.LinearQuadraticModel(GurobiSolver())

# Max           z
# s.t.  x + y     >= 0
#       x + y     >= 0
#      2x - y     >= 0
#       x     - z >= 0
#               z >= -10
# x ∈ [0, 1], y ∈ [0, ∞), z \in (-∞, 1]
MathProgBase.loadproblem!(m, [1 1 0; 1 1 0; 2 -1 0; 1 0 -1; 0 0 1], [0, 0, -Inf],[1,Inf,1], [0,0,1], [0,0,0,0,-10], [Inf,Inf,Inf,Inf,Inf], :Max)

MathProgBase.optimize!(m)

MathProgBase.addconstr!(m, [1, 2, 3], [-1.0000019, 1., -1.], 0., Inf)

MathProgBase.optimize!(m)

obj = MathProgBase.getobj(m)
MathProgBase.setobj!(m, obj)

MathProgBase.optimize!(m)
