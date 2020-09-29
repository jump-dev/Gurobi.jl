using Gurobi
using Test

# When adding new tests, be mindful about creating new environments. Either
# re-use an existing environment in the module, or name the test function
# `test_MULTI_ENV_xxx` to trap the specific Gurobi error indicating that an
# environment could not be created.

const GRB_ENV = Gurobi.Env()

function runtests(mod)
    for name in names(mod; all = true)
        sname = "$(name)"
        if !startswith(sname, "test_")
            continue
        end
        @testset "$(name)" begin
            if startswith(sname, "test_MULTI_ENV")
                try
                    getfield(mod, name)()
                catch ex
                    if ex == ErrorException(
                        "Gurobi Error 10009: Failed to obtain a valid license"
                    )
                        @warn(
                            "Skipping a test because there was an issue " *
                            "creating multiple licenses. This is probably " *
                            "because you have a limited license."
                        )
                    else
                        rethrow(ex)
                    end
                end
            else
                getfield(mod, name)()
            end
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
