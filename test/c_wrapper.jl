using Gurobi: getcoeff

@testset "File Tests" begin
    @testset for file in filter(
            s -> endswith(s, ".jl"), readdir(joinpath(@__DIR__, "C")))
        evalfile(joinpath(@__DIR__, "C", file))
    end
end

@testset "Empty constraints (Issue #142)" begin
    @testset "No variables, no constraints" begin
        model = Gurobi.Model(Gurobi.Env(), "model")
        A = Gurobi.get_constrmatrix(model)
        @test size(A) == (0, 0)
    end
    @testset "One variable, no constraints" begin
        model = Gurobi.Model(Gurobi.Env(), "model")
        Gurobi.add_cvar!(model, 0.0)
        Gurobi.update_model!(model)
        A = Gurobi.get_constrmatrix(model)
        @test size(A) == (0, 1)
    end
end

@testset "changing coefficients" begin
    A = [1.0 2.0; 3.0 0.0]
    b = [0.1, 0.2]
    f = [0.0, 0.0]
    model = gurobi_model(Gurobi.Env(); f=f, A=A, b=b)

    @test getcoeff(model, 1, 1) == 1.0

    # Verify that we can pass any kind of `::Integer` to `getcoeff`
    @test getcoeff(model, 1, Int8(1)) == 1.0
    @test getcoeff(model, UInt32(1), UInt32(1)) == 1.0

    @test getcoeff(model, 1, 2) == 2.0
    @test getcoeff(model, 2, 1) == 3.0
    @test getcoeff(model, 2, 2) == 0.0

    chg_coeffs!(model, 2, 1, 1.5)
    # before updating the model, we still get the old coefficient value
    @test getcoeff(model, 2, 1) == 3.0
    # after updating the model, we see the new coefficient value
    update_model!(model)
    @test getcoeff(model, 2, 1) == 1.5

    # Verify that we are automatically converting the scalars to
    # Cint and Float64 as necessary
    chg_coeffs!(model, Int8(2), Int128(1), Float32(2.0))
    update_model!(model)
    @test getcoeff(model, 2, 1) == 2.0

    chg_coeffs!(model, [1, 2], [1, 1], [5.0, 6.5])
    update_model!(model)
    @test getcoeff(model, 1, 1) == 5.0
    @test getcoeff(model, 2, 1) == 6.5
    @test getcoeff(model, 1, 2) == 2.0
    @test getcoeff(model, 2, 2) == 0.0

    chg_coeffs!(model, Int8[1, 2], Int8[1, 1], Float32[1.0, 2.0])
    update_model!(model)
    @test getcoeff(model, 1, 1) == 1.0
    @test getcoeff(model, 2, 1) == 2.0
    @test getcoeff(model, 1, 2) == 2.0
    @test getcoeff(model, 2, 2) == 0.0
end

@testset "Unicode support for names" begin
    model = Gurobi.Model(Gurobi.Env(), "model")
    Gurobi.add_cvar!(model, 0.0)
    for str in [
        "αβ", "Aα", "αA", "x̂", "xx̂", "x̂x", "¹", "x¹", "¹x", "₂", "x₂", "₂x",
        "xß", "̂", "̄"
    ]
        Gurobi.set_strattr!(model, "ModelName", str)
        Gurobi.update_model!(model)
        @test str == Gurobi.get_strattr(model, "ModelName")

        Gurobi.set_strattrelement!(model, "VarName", 1, str)
        Gurobi.update_model!(model)
        @test str == Gurobi.get_strattrelement(model, "VarName", 1)
    end
end
