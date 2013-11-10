using Gurobi

include(joinpath(Pkg.dir("MathProgBase"),"test","linprog.jl"))
linprogtest(GurobiSolver())

include(joinpath(Pkg.dir("MathProgBase"),"test","linproginterface.jl"))
linprogsolvertest(GurobiSolver())

include(joinpath(Pkg.dir("MathProgBase"),"test","mixintprog.jl"))
mixintprogtest(GurobiSolver())
