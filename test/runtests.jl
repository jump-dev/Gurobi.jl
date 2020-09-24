using Gurobi
using Random
using Test

@testset "MathOptInterface Tests" begin
    @testset "$(file)" for file in readdir("MOI")
        include(joinpath("MOI", file))
    end
end

@testset "Deprecated functions" begin
    err = ErrorException(Gurobi._DEPRECATED_ERROR_MESSAGE)
    @test_throws err get_status()
    @test_throws err Gurobi.get_status_code()
end
