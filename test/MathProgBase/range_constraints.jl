using Gurobi, MathProgBase, Test

@testset "Range Constraints" begin

    # Test Range constraints pass
    m = MathProgBase.LinearQuadraticModel(GurobiSolver(OutputFlag=0))
    MathProgBase.loadproblem!(m, [1 1], [0, 0],[1,1], [1,1], [0], [1], :Max)
    MathProgBase.optimize!(m)

    @test MathProgBase.getsolution(m) == [0.0, 1.0, 0.0]
    @test MathProgBase.getobjval(m) == 1.0
end
