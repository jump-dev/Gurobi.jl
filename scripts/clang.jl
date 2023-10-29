# Copyright (c) 2015 Dahua Lin, Miles Lubin, Joey Huchette, Iain Dunning, and
#   contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

# TODO(odow):
#
# This script can be used to build the C interface to Gurobi. However, it
# requires you to manually set the path to the appropriate gurobi_c.h.

# Note that clang 0.12 requires julia 1.5

import Clang

const LIBGRB_HEADERS = [
    "/Library/gurobi1100/macos_universal2/include/gurobi_c.h",
]

const GRB_VERSION = "110"

const GEN_DIR = joinpath(dirname(@__DIR__), "src", "gen$(GRB_VERSION)")
if !isdir(GEN_DIR)
    mkdir(GEN_DIR)
end

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
        lines[i] = replace(line, "Cstring" => "Ptr{Cchar}")
    end
    open(filename, "w") do io
        print.(Ref(io), lines)
    end
end
manual_corrections()

rm(joinpath(GEN_DIR, "LibTemplate.jl"))
