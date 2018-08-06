using Gurobi

using Compat
using Compat.Test, Compat.SparseArrays, Compat.Random

@testset "C API" begin
    include("c_wrapper.jl")
end

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
    for file in mpb_tests
        evalfile(joinpath("MathProgBase", "$(file).jl"))
    end
end

@testset "MathOptInterface Tests" begin
    include("MOIWrapper.jl")
end
