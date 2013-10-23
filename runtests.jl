tests = ["lp_01a", 
         "lp_01b", 
         "lp_02", 
         "lp_03",
         "mip_01", 
         "qp_01",
         "qp_02",
         "qcqp_01",
         "mathprog"]

for t in tests
    fp = joinpath("test", "$(t).jl")
    println("running $(fp) ...")
    include(fp)
end         
