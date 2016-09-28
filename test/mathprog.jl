using Gurobi

include(joinpath(Pkg.dir("MathProgBase"),"test","linprog.jl"))
linprogtest(GurobiSolver(InfUnbdInfo=1))

s = GurobiSolver()
MathProgBase.setparameters!(s, Silent=true, TimeLimit=100.0)

include(joinpath(Pkg.dir("MathProgBase"),"test","linproginterface.jl"))
linprogsolvertest(s)

include(joinpath(Pkg.dir("MathProgBase"),"test","mixintprog.jl"))
mixintprogtest(GurobiSolver())


include(joinpath(Pkg.dir("MathProgBase"),"test","quadprog.jl"))
quadprogtest(GurobiSolver())
socptest(GurobiSolver())
qpdualtest(GurobiSolver(QCPDual=1))

include(joinpath(Pkg.dir("MathProgBase"),"test","conicinterface.jl"))
coniclineartest(GurobiSolver())
conicSOCtest(GurobiSolver())
conicSOCRotatedtest(GurobiSolver())
conicSOCINTtest(GurobiSolver())
