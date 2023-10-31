# Copyright (c) 2015 Dahua Lin, Miles Lubin, Joey Huchette, Iain Dunning, and
#   contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

# !!! note
#
#     Run this script as `julia --project=. gen.jl`
#
#     When updating, you must:
#
#      * modify the `gurobi_c` constant to point to the correct filename
#      * modify the `output_folder`

using Clang.Generators

const gurobi_c = "/Library/gurobi1003/macos_universal2/include/gurobi_c.h"
const output_folder = "gen100"

options = load_options(joinpath(@__DIR__, "generate.toml"))
options["general"]["output_file_path"] =
    joinpath(@__DIR__, "..", "src", output_folder, "libgrb_api.jl")
build!(create_context([gurobi_c], get_default_args(), options))
