# TODO(odow):
#
# This script can be used to build the C interface to Gurobi. However, it requires
# you to manually do the following steps first:
#
# 1) Copy gurobi_c.h from Gurobi into this /scripts directory

import Clang

const LIBGRB_HEADERS = [
    joinpath(@__DIR__, "gurobi_c.h"),
]

const GEN_DIR = joinpath(dirname(@__DIR__), "src", "gen")

wc = Clang.init(
    headers = LIBGRB_HEADERS,
    output_file = joinpath(GEN_DIR, "libgrb_api.jl"),
    common_file = joinpath(GEN_DIR, "libgrb_common.jl"),
    header_wrapped = (root, current) -> root == current,
    header_library = x -> "libgurobi",
    clang_diagnostics = true,
)

run(wc)

function manual_corrections()
    filename = joinpath(GEN_DIR, "libgrb_api.jl")
    lines = readlines(filename; keep = true)
    for (i, line) in enumerate(lines)
        if occursin("GRBget", line) || occursin("GRBXget", line)
            lines[i] = replace(line, "Ptr{Cstring}" => "Ref{Ptr{Cchar}}")
            lines[i] = replace(line, "Cstring" => "Ptr{Cchar}")
        end
    end
    open(filename, "w") do io
        print.(Ref(io), lines)
    end
end
manual_corrections()

rm(joinpath(dirname(@__DIR__), "src", "gen", "LibTemplate.jl"))
