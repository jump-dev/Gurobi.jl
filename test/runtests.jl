using Gurobi
using Test

# When adding new tests, be mindful about creating new environments. Either
# re-use an existing environment in the module, or name the test function
# `test_MULTI_ENV_xxx` to trap the specific Gurobi error indicating that an
# environment could not be created.
const GRB_ENV = Gurobi.Env()

@testset "MathOptInterface Tests" begin
    for file in readdir("MOI")
        include(joinpath("MOI", file))
    end
end

@testset "Deprecated functions" begin
    err = ErrorException(Gurobi._DEPRECATED_ERROR_MESSAGE)
    @test_throws err get_status()
    @test_throws err Gurobi.get_status_code()
end
