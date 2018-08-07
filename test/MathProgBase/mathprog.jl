using Gurobi
import MathProgBase

if VERSION >= v"0.7-"
    using Compat.Pkg
end

function mathprogbase_file(file::String)
    if VERSION >= v"0.7-"
        return joinpath(dirname(dirname(pathof(MathProgBase))), "test", file)
    else
        return joinpath(Pkg.dir("MathProgBase"), "test", file)
    end
end

include(mathprogbase_file("linprog.jl"))
linprogtest(GurobiSolver(OutputFlag=0, InfUnbdInfo=1))

s = GurobiSolver()
MathProgBase.setparameters!(s, Silent=true, TimeLimit=100.0)

include(mathprogbase_file("linproginterface.jl"))
linprogsolvertest(s)

include(mathprogbase_file("mixintprog.jl"))
mixintprogtest(GurobiSolver(OutputFlag=0))


include(mathprogbase_file("quadprog.jl"))
quadprogtest(GurobiSolver(OutputFlag=0))
socptest(GurobiSolver(OutputFlag=0))
qpdualtest(GurobiSolver(OutputFlag=0, QCPDual=1))

include(mathprogbase_file("conicinterface.jl"))
coniclineartest(GurobiSolver(OutputFlag=0))
# The following tests are not passing on Gurobi 8.0.0 due to known bug
# on infeasibility check of SOC models.
# https://github.com/JuliaOpt/Gurobi.jl/pull/123
if Gurobi.version < v"8.0.0"
    conicSOCtest(GurobiSolver(OutputFlag=0))
    conicSOCRotatedtest(GurobiSolver(OutputFlag=0))
end
conicSOCINTtest(GurobiSolver(OutputFlag=0))
