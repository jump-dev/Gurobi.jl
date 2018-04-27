using Gurobi, Base.Test, MathOptInterface, MathOptInterface.Test

const MOIT = MathOptInterface.Test

@testset "Linear tests" begin
    linconfig = linconfig = MOIT.TestConfig()
    solver = GurobiOptimizer(OutputFlag=0)
    MOIT.contlineartest(solver, linconfig, ["linear10","linear12","linear8a","linear8b","linear8c"])

    solver_nopresolve = GurobiOptimizer(OutputFlag=0, InfUnbdInfo=1)
    MOIT.contlineartest(solver_nopresolve, linconfig, ["linear10","linear12","linear8a"])

    linconfig_nocertificate = MOIT.TestConfig(infeas_certificates=false)
    MOIT.linear12test(solver, linconfig_nocertificate)
    MOIT.linear8atest(solver, linconfig_nocertificate)

    # 10 is ranged
end

@testset "Quadratic tests" begin
    quadconfig = MOIT.TestConfig(atol=1e-4, rtol=1e-4, duals=false, query=false)
    solver = GurobiOptimizer(OutputFlag=0)
    MOIT.contquadratictest(solver, quadconfig)
end

@testset "Linear Conic tests" begin
    linconfig = MOIT.TestConfig()
    solver = GurobiOptimizer(OutputFlag=0)
    MOIT.lintest(solver, linconfig, ["lin3","lin4"])

    solver_nopresolve = GurobiOptimizer(OutputFlag=0, InfUnbdInfo=1)
    MOIT.lintest(solver_nopresolve, linconfig)
end

@testset "Integer Linear tests" begin
    intconfig = MOIT.TestConfig()
    solver = GurobiOptimizer(OutputFlag=0)
    MOIT.intlineartest(solver, intconfig, ["int3"])

    # 3 is ranged
end
@testset "ModelLike tests" begin
    intconfig = MOIT.TestConfig()
    solver = GurobiOptimizer()
    MOIT.validtest(solver)
    MOIT.emptytest(solver)
    solver2 = GurobiOptimizer()
    MOIT.copytest(solver,solver2)
end
