if haskey(ENV, "GITHUB_ACTIONS")
    # We're being run as part of a Github action. The most likely case is that
    # this is the auto-merge action as part of the General registry.
    # For now, we're going to silently skip the tests.
    @info("Detected a Github action. Skipping tests.")
    exit(0)
end

using Gurobi, Test, Random

@testset "C API" begin
    include("c_wrapper.jl")
end

@testset "MathProgBase Tests" begin
    @testset for file in ["env", "large_coefficients", "range_constraints",]
        evalfile(joinpath("MathProgBase", "$(file).jl"))
    end
    @testset "MathProgJuMP" begin
        evalfile(joinpath("MathProgBase", "mathprog.jl"))
    end
end

@testset "MathOptInterface Tests" begin
    @testset "$(file)" for file in readdir("MOI")
        include(joinpath("MOI", file))
    end
end
