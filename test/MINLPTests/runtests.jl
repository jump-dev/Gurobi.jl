# Copyright (c) 2015 Dahua Lin, Miles Lubin, Joey Huchette, Iain Dunning, and
#   contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

using Gurobi
using MINLPTests
using Test

@testset "MINLPTests" begin
    tol = 1e-2
    optimizer = MINLPTests.JuMP.optimizer_with_attributes(
        Gurobi.Optimizer,
        "TimeLimit" => 60.0,
    )
    MINLPTests.test_nlp_cvx_expr(
        optimizer;
        termination_target = MINLPTests.TERMINATION_TARGET_GLOBAL,
        primal_target = MINLPTests.PRIMAL_TARGET_GLOBAL,
        objective_tol = 1e-6,
        primal_tol = tol,
        dual_tol = tol,
    )
    MINLPTests.test_nlp_expr(
        optimizer;
        termination_target = MINLPTests.TERMINATION_TARGET_GLOBAL,
        primal_target = MINLPTests.PRIMAL_TARGET_GLOBAL,
        objective_tol = 1e-6,
        primal_tol = tol,
        dual_tol = tol,
    )
    MINLPTests.test_nlp_mi_expr(
        optimizer;
        termination_target = MINLPTests.TERMINATION_TARGET_GLOBAL,
        primal_target = MINLPTests.PRIMAL_TARGET_GLOBAL,
        objective_tol = 1e-6,
        primal_tol = tol,
        dual_tol = tol,
    )
end