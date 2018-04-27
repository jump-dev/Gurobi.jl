using Base.Test

const mpb_tests = [
    "lp_01a",
    "lp_01b",
    "lp_02",
    "lp_03",
    "lp_04",
    "mip_01",
    "qp_01",
    "qp_02",
    "qcqp_01",
    "mathprog",
    "test_grb_attrs",
    "env",
    "range_constraints",
    "test_get_strarray",
    "large_coefficients",
    "multiobj",
    "test_read"
]

@testset "MathProgBase Tests" begin
    for t in mpb_tests
        fp = "$(t).jl"
        println("running $(fp) ...")
        evalfile(joinpath("MathProgBase", fp))
    end
end

const moi_tests = [
    "mathoptinterface_tests"
]

@testset "MathOptInterface Tests" begin
    for t in moi_tests
        fp = "$(t).jl"
        println("running $(fp) ...")
        evalfile(joinpath("MathOptInterface", fp))
    end
end
