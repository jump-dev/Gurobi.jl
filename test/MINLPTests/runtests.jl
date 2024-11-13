# Copyright (c) 2015 Dahua Lin, Miles Lubin, Joey Huchette, Iain Dunning, and
#   contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

using Gurobi
using MINLPTests
using Tests

@testset "MINLPTests" begin
    MINLPTests.test_nlp_cvx_expr(
        Gurobi.Optimizer;
        termination_target = MINLPTests.TERMINATION_TARGET_GLOBAL,
        primal_target = MINLPTests.PRIMAL_TARGET_GLOBAL,
        objective_tol = 1e-6,
        primal_tol = 1e-6,
        dual_tol = 1e-6,
    )

    MINLPTests.test_nlp_expr(
        Gurobi.Optimizer;
        termination_target = MINLPTests.TERMINATION_TARGET_GLOBAL,
        primal_target = MINLPTests.PRIMAL_TARGET_GLOBAL,
        objective_tol = 1e-6,
        primal_tol = 1e-6,
        dual_tol = 1e-6,
    )

    MINLPTests.test_nlp_mi_expr(
        Gurobi.Optimizer;
        termination_target = MINLPTests.TERMINATION_TARGET_GLOBAL,
        primal_target = MINLPTests.PRIMAL_TARGET_GLOBAL,
        objective_tol = 1e-6,
        primal_tol = 1e-6,
        dual_tol = 1e-6,
    )
end