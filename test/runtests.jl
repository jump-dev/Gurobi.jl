using Gurobi
using Test

function runtests(mod)
    for name in names(mod; all = true)
        if !startswith("$(name)", "test_")
            continue
        end
        @testset "$(name)" begin
            getfield(mod, name)()
        end
    end
end

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
