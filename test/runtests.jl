# Copyright (c) 2015 Dahua Lin, Miles Lubin, Joey Huchette, Iain Dunning, and
#   contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

using Gurobi
using Test

# When adding new tests, be mindful about creating new environments. Either
# re-use an existing environment in the module, or name the test function
# `test_MULTI_ENV_xxx` to trap the specific Gurobi error indicating that an
# environment could not be created.
const GRB_ENV = Gurobi.Env(output_flag = 0)

@testset "MathOptInterface Tests" begin
    for file in readdir("MOI")
        include(joinpath("MOI", file))
    end
end
