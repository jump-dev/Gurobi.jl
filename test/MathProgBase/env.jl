using Gurobi, MathProgBase, Test

@testset "Env" begin
    env = Gurobi.Env()

    # Test supplying environment manually

    solver = GurobiSolver(env)
    # Check that the env for each model is the same
    m1 = MathProgBase.LinearQuadraticModel(solver)
    @test m1.inner.env === env
    m2 = MathProgBase.LinearQuadraticModel(solver)
    @test m2.inner.env === env
    # Check that finalizer doesn't touch env when manually provided
    finalize(m1.inner)
    @test Gurobi.is_valid(env)

    # Test creating environment automatically

    solver = GurobiSolver()
    # Check that the env for each model is different
    m3 = MathProgBase.LinearQuadraticModel(solver)
    @test m3.inner.env !== env
    m4 = MathProgBase.LinearQuadraticModel(solver)
    @test m4.inner.env !== env
    @test m3.inner.env !== m4.inner.env
    # Check that env is finalized with model when not supplied manually
    finalize(m3.inner)
    @test !Gurobi.is_valid(m3.inner.env)
end
