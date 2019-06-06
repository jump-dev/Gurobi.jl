using Gurobi, MathProgBase, Test

@testset "Large Coefficients" begin

    @testset "Large Objective Coefficients" begin
        # Min   1.1e100x
        #       x >= 1
        m = MathProgBase.LinearQuadraticModel(GurobiSolver(OutputFlag=0))
        MathProgBase.loadproblem!(m, Array{Float64}(undef, 0, 1), [1],[Inf], [1.1e100], Float64[], Float64[], :Min)
        MathProgBase.optimize!(m)
        @test MathProgBase.getsolution(m) == [1.0]
        @test MathProgBase.getobjval(m) == 1e100

        # Max   -1.1e100x
        #       x >= 1
        m = MathProgBase.LinearQuadraticModel(GurobiSolver(OutputFlag=0))
        MathProgBase.loadproblem!(m, Array{Float64}(undef, 0, 1), [1],[Inf], [-1.1e100], Float64[], Float64[], :Max)
        MathProgBase.optimize!(m)
        @test MathProgBase.getsolution(m) == [1.0]
        @test MathProgBase.getobjval(m) == -1e100
    end

    @testset "Large Variable Bounds" begin
        # Min   x
        #       x >= 1.1e30
        m = MathProgBase.LinearQuadraticModel(GurobiSolver(OutputFlag=0))
        MathProgBase.loadproblem!(m, Array{Float64}(undef, 0, 1), [1.1e30],[Inf], [1], Float64[], Float64[], :Min)
        MathProgBase.optimize!(m)
        @test MathProgBase.status(m) == :Infeasible

        # Max   x
        #       x <= -1.1e30
        m = MathProgBase.LinearQuadraticModel(GurobiSolver(OutputFlag=0))
        MathProgBase.loadproblem!(m, Array{Float64}(undef, 0, 1), [-Inf],[-1.1e30], [1], Float64[], Float64[], :Max)
        MathProgBase.optimize!(m)
        @test MathProgBase.status(m) == :Infeasible

        # Min   x
        #       x >= -1.1e30
        m = MathProgBase.LinearQuadraticModel(GurobiSolver(OutputFlag=0))
        MathProgBase.loadproblem!(m, Array{Float64}(undef, 0, 1), [-1.1e30],[Inf], [1], Float64[], Float64[], :Min)
        MathProgBase.optimize!(m)
        @test MathProgBase.status(m) == :Unbounded

        # Max   x
        #       x <= 1.1e30
        m = MathProgBase.LinearQuadraticModel(GurobiSolver(OutputFlag=0))
        MathProgBase.loadproblem!(m, Array{Float64}(undef, 0, 1), [-Inf],[1.1e30], [1], Float64[], Float64[], :Max)
        MathProgBase.optimize!(m)
        @test MathProgBase.status(m) == :Unbounded

    end
end
