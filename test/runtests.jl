tests = [
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

for t in tests
    fp = "$(t).jl"
    println("running $(fp) ...")
    evalfile(fp)
end
