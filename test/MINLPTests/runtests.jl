# Copyright (c) 2015 Dahua Lin, Miles Lubin, Joey Huchette, Iain Dunning, and
#   contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

using Gurobi
using MINLPTests
using Test

@testset "MINLPTests" begin
    optimizer = MINLPTests.JuMP.optimizer_with_attributes(
        Gurobi.Optimizer,
        "TimeLimit" => 60.0,
    )
    MINLPTests.test_nlp_cvx_expr(
        optimizer;
        termination_target = MINLPTests.TERMINATION_TARGET_GLOBAL,
        primal_target = MINLPTests.PRIMAL_TARGET_GLOBAL,
        objective_tol = 1e-6,
        primal_tol = 1e-2,
        dual_tol = NaN,
        exclude = [
            "501_011",  # TIME_LIMIT
        ],
    )
    MINLPTests.test_nlp_expr(
        optimizer;
        termination_target = MINLPTests.TERMINATION_TARGET_GLOBAL,
        primal_target = MINLPTests.PRIMAL_TARGET_GLOBAL,
        objective_tol = 1e-6,
        primal_tol = 1e-2,
        dual_tol = NaN,
        exclude = [
            "001_010",  # TIME_LIMIT
            "004_010",  # abs
            "004_011",  # abs
            "005_010",  # inv
            "006_010",  # UserDefinedFunction
            "009_010",  # min
            "009_011",  # max
        ],
    )
    MINLPTests.test_nlp_mi_expr(
        optimizer;
        termination_target = MINLPTests.TERMINATION_TARGET_GLOBAL,
        primal_target = MINLPTests.PRIMAL_TARGET_GLOBAL,
        objective_tol = 1e-6,
        primal_tol = 1e-2,
        dual_tol = NaN,
        exclude = [
            "001_010",  # TIME_LIMIT
            "004_010",  # abs
            "004_011",  # abs
            "005_010",  # inv
            "006_010",  # UserDefinedFunction
        ],
    )
end
