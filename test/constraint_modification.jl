using Gurobi
using Base.Test
using Gurobi: getcoeff


@testset "changing coefficients" begin
    A = spzeros(2, 2)
    A[1, 1] = 1.0
    A[1, 2] = 2.0
    A[2, 1] = 3.0
    b = [0.1, 0.2]
    f = [0.0, 0.0]
    model = gurobi_model(Gurobi.Env(); f=f, A=A, b=b)

    @test getcoeff(model, 1, 1) == 1.0
    @test getcoeff(model, 1, 2) == 2.0
    @test getcoeff(model, 2, 1) == 3.0
    @test getcoeff(model, 2, 2) == 0.0

    chg_coeffs!(model, 2, 1, 1.5)
    # before updating the model, we still get the old coefficient value
    @test getcoeff(model, 2, 1) == 3.0
    # after updating the model, we see the new coefficient value
    update_model!(model)
    @test getcoeff(model, 2, 1) == 1.5

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
