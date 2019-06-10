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
    include("MOI_wrapper.jl")
end
