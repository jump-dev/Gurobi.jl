using Gurobi
using Random
using Test

@testset "MathOptInterface Tests" begin
    @testset "$(file)" for file in readdir("MOI")
        include(joinpath("MOI", file))
    end
end
