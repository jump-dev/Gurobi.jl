# Copyright (c) 2015 Dahua Lin, Miles Lubin, Joey Huchette, Iain Dunning, and
#   contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

# !format: off

# ============================ start of prologue.jl ============================
# These constants are handled explicitly in the prologue.jl to avoid them being
# automatically parsed as `const GRB_LESS_EQUAL = Cchar('<')`. There's probably
# a better way to handle this, but it works for now, and they won't be changing
# in future releases.
const GRB_LESS_EQUAL = '<'
const GRB_GREATER_EQUAL = '>'
const GRB_EQUAL = '='
const GRB_CONTINUOUS = 'C'
const GRB_BINARY = 'B'
const GRB_INTEGER = 'I'
const GRB_SEMICONT = 'S'
const GRB_SEMIINT = 'N'
# ============================= end of prologue.jl =============================
