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

const _GRBmodel = Cvoid

const GRBmodel = _GRBmodel

const _GRBenv = Cvoid

const GRBenv = _GRBenv

function GRBemptyenvadvinternal(
    envP,
    apitype,
    major,
    minor,
    tech,
    malloccb,
    calloccb,
    realloccb,
    freecb,
    threadcreatecb,
    threadjoincb,
    syscbusrdata,
)
    return ccall(
        (:GRBemptyenvadvinternal, libgurobi),
        Cint,
        (
            Ptr{Ptr{GRBenv}},
            Cint,
            Cint,
            Cint,
            Cint,
            Ptr{Cvoid},
            Ptr{Cvoid},
            Ptr{Cvoid},
            Ptr{Cvoid},
            Ptr{Cvoid},
            Ptr{Cvoid},
            Ptr{Cvoid},
        ),
        envP,
        apitype,
        major,
        minor,
        tech,
        malloccb,
        calloccb,
        realloccb,
        freecb,
        threadcreatecb,
        threadjoincb,
        syscbusrdata,
    )
end

function GRBloadenvinternal(envP, logfilename, major, minor, tech)
    return ccall(
        (:GRBloadenvinternal, libgurobi),
        Cint,
        (Ptr{Ptr{GRBenv}}, Ptr{Cchar}, Cint, Cint, Cint),
        envP,
        logfilename,
        major,
        minor,
        tech,
    )
end

function GRBemptyenvinternal(envP, major, minor, tech)
    return ccall(
        (:GRBemptyenvinternal, libgurobi),
        Cint,
        (Ptr{Ptr{GRBenv}}, Cint, Cint, Cint),
        envP,
        major,
        minor,
        tech,
    )
end

const _GRBbatch = Cvoid

const GRBbatch = _GRBbatch

function GRBgetattrinfo(model, attrname, datatypeP, attrtypeP, settableP)
    return ccall(
        (:GRBgetattrinfo, libgurobi),
        Cint,
        (Ptr{GRBmodel}, Ptr{Cchar}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}),
        model,
        attrname,
        datatypeP,
        attrtypeP,
        settableP,
    )
end

function GRBisattravailable(model, attrname)
    return ccall(
        (:GRBisattravailable, libgurobi),
        Cint,
        (Ptr{GRBmodel}, Ptr{Cchar}),
        model,
        attrname,
    )
end

function GRBgetintattr(model, attrname, valueP)
    return ccall(
        (:GRBgetintattr, libgurobi),
        Cint,
        (Ptr{GRBmodel}, Ptr{Cchar}, Ptr{Cint}),
        model,
        attrname,
        valueP,
    )
end

function GRBsetintattr(model, attrname, newvalue)
    return ccall(
        (:GRBsetintattr, libgurobi),
        Cint,
        (Ptr{GRBmodel}, Ptr{Cchar}, Cint),
        model,
        attrname,
        newvalue,
    )
end

function GRBgetintattrelement(model, attrname, element, valueP)
    return ccall(
        (:GRBgetintattrelement, libgurobi),
        Cint,
        (Ptr{GRBmodel}, Ptr{Cchar}, Cint, Ptr{Cint}),
        model,
        attrname,
        element,
        valueP,
    )
end

function GRBsetintattrelement(model, attrname, element, newvalue)
    return ccall(
        (:GRBsetintattrelement, libgurobi),
        Cint,
        (Ptr{GRBmodel}, Ptr{Cchar}, Cint, Cint),
        model,
        attrname,
        element,
        newvalue,
    )
end

function GRBgetintattrarray(model, attrname, first, len, values)
    return ccall(
        (:GRBgetintattrarray, libgurobi),
        Cint,
        (Ptr{GRBmodel}, Ptr{Cchar}, Cint, Cint, Ptr{Cint}),
        model,
        attrname,
        first,
        len,
        values,
    )
end

function GRBsetintattrarray(model, attrname, first, len, newvalues)
    return ccall(
        (:GRBsetintattrarray, libgurobi),
        Cint,
        (Ptr{GRBmodel}, Ptr{Cchar}, Cint, Cint, Ptr{Cint}),
        model,
        attrname,
        first,
        len,
        newvalues,
    )
end

function GRBgetintattrlist(model, attrname, len, ind, values)
    return ccall(
        (:GRBgetintattrlist, libgurobi),
        Cint,
        (Ptr{GRBmodel}, Ptr{Cchar}, Cint, Ptr{Cint}, Ptr{Cint}),
        model,
        attrname,
        len,
        ind,
        values,
    )
end

function GRBsetintattrlist(model, attrname, len, ind, newvalues)
    return ccall(
        (:GRBsetintattrlist, libgurobi),
        Cint,
        (Ptr{GRBmodel}, Ptr{Cchar}, Cint, Ptr{Cint}, Ptr{Cint}),
        model,
        attrname,
        len,
        ind,
        newvalues,
    )
end

function GRBgetcharattrelement(model, attrname, element, valueP)
    return ccall(
        (:GRBgetcharattrelement, libgurobi),
        Cint,
        (Ptr{GRBmodel}, Ptr{Cchar}, Cint, Ptr{Cchar}),
        model,
        attrname,
        element,
        valueP,
    )
end

function GRBsetcharattrelement(model, attrname, element, newvalue)
    return ccall(
        (:GRBsetcharattrelement, libgurobi),
        Cint,
        (Ptr{GRBmodel}, Ptr{Cchar}, Cint, Cchar),
        model,
        attrname,
        element,
        newvalue,
    )
end

function GRBgetcharattrarray(model, attrname, first, len, values)
    return ccall(
        (:GRBgetcharattrarray, libgurobi),
        Cint,
        (Ptr{GRBmodel}, Ptr{Cchar}, Cint, Cint, Ptr{Cchar}),
        model,
        attrname,
        first,
        len,
        values,
    )
end

function GRBsetcharattrarray(model, attrname, first, len, newvalues)
    return ccall(
        (:GRBsetcharattrarray, libgurobi),
        Cint,
        (Ptr{GRBmodel}, Ptr{Cchar}, Cint, Cint, Ptr{Cchar}),
        model,
        attrname,
        first,
        len,
        newvalues,
    )
end

function GRBgetcharattrlist(model, attrname, len, ind, values)
    return ccall(
        (:GRBgetcharattrlist, libgurobi),
        Cint,
        (Ptr{GRBmodel}, Ptr{Cchar}, Cint, Ptr{Cint}, Ptr{Cchar}),
        model,
        attrname,
        len,
        ind,
        values,
    )
end

function GRBsetcharattrlist(model, attrname, len, ind, newvalues)
    return ccall(
        (:GRBsetcharattrlist, libgurobi),
        Cint,
        (Ptr{GRBmodel}, Ptr{Cchar}, Cint, Ptr{Cint}, Ptr{Cchar}),
        model,
        attrname,
        len,
        ind,
        newvalues,
    )
end

function GRBgetdblattr(model, attrname, valueP)
    return ccall(
        (:GRBgetdblattr, libgurobi),
        Cint,
        (Ptr{GRBmodel}, Ptr{Cchar}, Ptr{Cdouble}),
        model,
        attrname,
        valueP,
    )
end

function GRBsetdblattr(model, attrname, newvalue)
    return ccall(
        (:GRBsetdblattr, libgurobi),
        Cint,
        (Ptr{GRBmodel}, Ptr{Cchar}, Cdouble),
        model,
        attrname,
        newvalue,
    )
end

function GRBgetdblattrelement(model, attrname, element, valueP)
    return ccall(
        (:GRBgetdblattrelement, libgurobi),
        Cint,
        (Ptr{GRBmodel}, Ptr{Cchar}, Cint, Ptr{Cdouble}),
        model,
        attrname,
        element,
        valueP,
    )
end

function GRBsetdblattrelement(model, attrname, element, newvalue)
    return ccall(
        (:GRBsetdblattrelement, libgurobi),
        Cint,
        (Ptr{GRBmodel}, Ptr{Cchar}, Cint, Cdouble),
        model,
        attrname,
        element,
        newvalue,
    )
end

function GRBgetdblattrarray(model, attrname, first, len, values)
    return ccall(
        (:GRBgetdblattrarray, libgurobi),
        Cint,
        (Ptr{GRBmodel}, Ptr{Cchar}, Cint, Cint, Ptr{Cdouble}),
        model,
        attrname,
        first,
        len,
        values,
    )
end

function GRBsetdblattrarray(model, attrname, first, len, newvalues)
    return ccall(
        (:GRBsetdblattrarray, libgurobi),
        Cint,
        (Ptr{GRBmodel}, Ptr{Cchar}, Cint, Cint, Ptr{Cdouble}),
        model,
        attrname,
        first,
        len,
        newvalues,
    )
end

function GRBgetdblattrlist(model, attrname, len, ind, values)
    return ccall(
        (:GRBgetdblattrlist, libgurobi),
        Cint,
        (Ptr{GRBmodel}, Ptr{Cchar}, Cint, Ptr{Cint}, Ptr{Cdouble}),
        model,
        attrname,
        len,
        ind,
        values,
    )
end

function GRBsetdblattrlist(model, attrname, len, ind, newvalues)
    return ccall(
        (:GRBsetdblattrlist, libgurobi),
        Cint,
        (Ptr{GRBmodel}, Ptr{Cchar}, Cint, Ptr{Cint}, Ptr{Cdouble}),
        model,
        attrname,
        len,
        ind,
        newvalues,
    )
end

function GRBgetstrattr(model, attrname, valueP)
    return ccall(
        (:GRBgetstrattr, libgurobi),
        Cint,
        (Ptr{GRBmodel}, Ptr{Cchar}, Ptr{Ptr{Cchar}}),
        model,
        attrname,
        valueP,
    )
end

function GRBsetstrattr(model, attrname, newvalue)
    return ccall(
        (:GRBsetstrattr, libgurobi),
        Cint,
        (Ptr{GRBmodel}, Ptr{Cchar}, Ptr{Cchar}),
        model,
        attrname,
        newvalue,
    )
end

function GRBgetstrattrelement(model, attrname, element, valueP)
    return ccall(
        (:GRBgetstrattrelement, libgurobi),
        Cint,
        (Ptr{GRBmodel}, Ptr{Cchar}, Cint, Ptr{Ptr{Cchar}}),
        model,
        attrname,
        element,
        valueP,
    )
end

function GRBsetstrattrelement(model, attrname, element, newvalue)
    return ccall(
        (:GRBsetstrattrelement, libgurobi),
        Cint,
        (Ptr{GRBmodel}, Ptr{Cchar}, Cint, Ptr{Cchar}),
        model,
        attrname,
        element,
        newvalue,
    )
end

function GRBgetstrattrarray(model, attrname, first, len, values)
    return ccall(
        (:GRBgetstrattrarray, libgurobi),
        Cint,
        (Ptr{GRBmodel}, Ptr{Cchar}, Cint, Cint, Ptr{Ptr{Cchar}}),
        model,
        attrname,
        first,
        len,
        values,
    )
end

function GRBsetstrattrarray(model, attrname, first, len, newvalues)
    return ccall(
        (:GRBsetstrattrarray, libgurobi),
        Cint,
        (Ptr{GRBmodel}, Ptr{Cchar}, Cint, Cint, Ptr{Ptr{Cchar}}),
        model,
        attrname,
        first,
        len,
        newvalues,
    )
end

function GRBgetstrattrlist(model, attrname, len, ind, values)
    return ccall(
        (:GRBgetstrattrlist, libgurobi),
        Cint,
        (Ptr{GRBmodel}, Ptr{Cchar}, Cint, Ptr{Cint}, Ptr{Ptr{Cchar}}),
        model,
        attrname,
        len,
        ind,
        values,
    )
end

function GRBsetstrattrlist(model, attrname, len, ind, newvalues)
    return ccall(
        (:GRBsetstrattrlist, libgurobi),
        Cint,
        (Ptr{GRBmodel}, Ptr{Cchar}, Cint, Ptr{Cint}, Ptr{Ptr{Cchar}}),
        model,
        attrname,
        len,
        ind,
        newvalues,
    )
end

function GRBsetcallbackfunc(model, cb, usrdata)
    return ccall(
        (:GRBsetcallbackfunc, libgurobi),
        Cint,
        (Ptr{GRBmodel}, Ptr{Cvoid}, Ptr{Cvoid}),
        model,
        cb,
        usrdata,
    )
end

function GRBgetcallbackfuncenv(env, cbP)
    return ccall(
        (:GRBgetcallbackfuncenv, libgurobi),
        Cint,
        (Ptr{GRBenv}, Ptr{Ptr{Cvoid}}),
        env,
        cbP,
    )
end

function GRBsetcallbackfuncenv(env, cb, usrdata)
    return ccall(
        (:GRBsetcallbackfuncenv, libgurobi),
        Cint,
        (Ptr{GRBenv}, Ptr{Cvoid}, Ptr{Cvoid}),
        env,
        cb,
        usrdata,
    )
end

function GRBgetcallbackfunc(model, cbP)
    return ccall(
        (:GRBgetcallbackfunc, libgurobi),
        Cint,
        (Ptr{GRBmodel}, Ptr{Ptr{Cvoid}}),
        model,
        cbP,
    )
end

function GRBsetlogcallbackfunc(model, cb, logdata)
    return ccall(
        (:GRBsetlogcallbackfunc, libgurobi),
        Cint,
        (Ptr{GRBmodel}, Ptr{Cvoid}, Ptr{Cvoid}),
        model,
        cb,
        logdata,
    )
end

function GRBsetlogcallbackfuncenv(env, cb, logdata)
    return ccall(
        (:GRBsetlogcallbackfuncenv, libgurobi),
        Cint,
        (Ptr{GRBenv}, Ptr{Cvoid}, Ptr{Cvoid}),
        env,
        cb,
        logdata,
    )
end

function GRBgetlogcallbackfuncenv(env, cbP, logdataP)
    return ccall(
        (:GRBgetlogcallbackfuncenv, libgurobi),
        Cint,
        (Ptr{GRBenv}, Ptr{Ptr{Cvoid}}, Ptr{Ptr{Cvoid}}),
        env,
        cbP,
        logdataP,
    )
end

function GRBcbproceed(cbdata)
    return ccall((:GRBcbproceed, libgurobi), Cint, (Ptr{Cvoid},), cbdata)
end

function GRBcbget(cbdata, where, what, resultP)
    return ccall(
        (:GRBcbget, libgurobi),
        Cint,
        (Ptr{Cvoid}, Cint, Cint, Ptr{Cvoid}),
        cbdata,
        where,
        what,
        resultP,
    )
end

function GRBcbsetintparam(cbdata, paramname, newvalue)
    return ccall(
        (:GRBcbsetintparam, libgurobi),
        Cint,
        (Ptr{Cvoid}, Ptr{Cchar}, Cint),
        cbdata,
        paramname,
        newvalue,
    )
end

function GRBcbsetdblparam(cbdata, paramname, newvalue)
    return ccall(
        (:GRBcbsetdblparam, libgurobi),
        Cint,
        (Ptr{Cvoid}, Ptr{Cchar}, Cdouble),
        cbdata,
        paramname,
        newvalue,
    )
end

function GRBcbsetstrparam(cbdata, paramname, newvalue)
    return ccall(
        (:GRBcbsetstrparam, libgurobi),
        Cint,
        (Ptr{Cvoid}, Ptr{Cchar}, Ptr{Cchar}),
        cbdata,
        paramname,
        newvalue,
    )
end

function GRBcbsetparam(cbdata, paramname, newvalue)
    return ccall(
        (:GRBcbsetparam, libgurobi),
        Cint,
        (Ptr{Cvoid}, Ptr{Cchar}, Ptr{Cchar}),
        cbdata,
        paramname,
        newvalue,
    )
end

function GRBcbsolution(cbdata, solution, objvalP)
    return ccall(
        (:GRBcbsolution, libgurobi),
        Cint,
        (Ptr{Cvoid}, Ptr{Cdouble}, Ptr{Cdouble}),
        cbdata,
        solution,
        objvalP,
    )
end

function GRBcbcut(cbdata, cutlen, cutind, cutval, cutsense, cutrhs)
    return ccall(
        (:GRBcbcut, libgurobi),
        Cint,
        (Ptr{Cvoid}, Cint, Ptr{Cint}, Ptr{Cdouble}, Cchar, Cdouble),
        cbdata,
        cutlen,
        cutind,
        cutval,
        cutsense,
        cutrhs,
    )
end

function GRBcblazy(cbdata, lazylen, lazyind, lazyval, lazysense, lazyrhs)
    return ccall(
        (:GRBcblazy, libgurobi),
        Cint,
        (Ptr{Cvoid}, Cint, Ptr{Cint}, Ptr{Cdouble}, Cchar, Cdouble),
        cbdata,
        lazylen,
        lazyind,
        lazyval,
        lazysense,
        lazyrhs,
    )
end

function GRBgetcoeff(model, constr, var, valP)
    return ccall(
        (:GRBgetcoeff, libgurobi),
        Cint,
        (Ptr{GRBmodel}, Cint, Cint, Ptr{Cdouble}),
        model,
        constr,
        var,
        valP,
    )
end

function GRBgetconstrs(model, numnzP, cbeg, cind, cval, start, len)
    return ccall(
        (:GRBgetconstrs, libgurobi),
        Cint,
        (
            Ptr{GRBmodel},
            Ptr{Cint},
            Ptr{Cint},
            Ptr{Cint},
            Ptr{Cdouble},
            Cint,
            Cint,
        ),
        model,
        numnzP,
        cbeg,
        cind,
        cval,
        start,
        len,
    )
end

function GRBXgetconstrs(model, numnzP, cbeg, cind, cval, start, len)
    return ccall(
        (:GRBXgetconstrs, libgurobi),
        Cint,
        (
            Ptr{GRBmodel},
            Ptr{Csize_t},
            Ptr{Csize_t},
            Ptr{Cint},
            Ptr{Cdouble},
            Cint,
            Cint,
        ),
        model,
        numnzP,
        cbeg,
        cind,
        cval,
        start,
        len,
    )
end

function GRBgetvars(model, numnzP, vbeg, vind, vval, start, len)
    return ccall(
        (:GRBgetvars, libgurobi),
        Cint,
        (
            Ptr{GRBmodel},
            Ptr{Cint},
            Ptr{Cint},
            Ptr{Cint},
            Ptr{Cdouble},
            Cint,
            Cint,
        ),
        model,
        numnzP,
        vbeg,
        vind,
        vval,
        start,
        len,
    )
end

function GRBXgetvars(model, numnzP, vbeg, vind, vval, start, len)
    return ccall(
        (:GRBXgetvars, libgurobi),
        Cint,
        (
            Ptr{GRBmodel},
            Ptr{Csize_t},
            Ptr{Csize_t},
            Ptr{Cint},
            Ptr{Cdouble},
            Cint,
            Cint,
        ),
        model,
        numnzP,
        vbeg,
        vind,
        vval,
        start,
        len,
    )
end

function GRBgetsos(model, nummembersP, sostype, beg, ind, weight, start, len)
    return ccall(
        (:GRBgetsos, libgurobi),
        Cint,
        (
            Ptr{GRBmodel},
            Ptr{Cint},
            Ptr{Cint},
            Ptr{Cint},
            Ptr{Cint},
            Ptr{Cdouble},
            Cint,
            Cint,
        ),
        model,
        nummembersP,
        sostype,
        beg,
        ind,
        weight,
        start,
        len,
    )
end

function GRBgetgenconstrMax(model, genconstr, resvarP, nvarsP, vars, constantP)
    return ccall(
        (:GRBgetgenconstrMax, libgurobi),
        Cint,
        (Ptr{GRBmodel}, Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}),
        model,
        genconstr,
        resvarP,
        nvarsP,
        vars,
        constantP,
    )
end

function GRBgetgenconstrMin(model, genconstr, resvarP, nvarsP, vars, constantP)
    return ccall(
        (:GRBgetgenconstrMin, libgurobi),
        Cint,
        (Ptr{GRBmodel}, Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}),
        model,
        genconstr,
        resvarP,
        nvarsP,
        vars,
        constantP,
    )
end

function GRBgetgenconstrAbs(model, genconstr, resvarP, argvarP)
    return ccall(
        (:GRBgetgenconstrAbs, libgurobi),
        Cint,
        (Ptr{GRBmodel}, Cint, Ptr{Cint}, Ptr{Cint}),
        model,
        genconstr,
        resvarP,
        argvarP,
    )
end

function GRBgetgenconstrAnd(model, genconstr, resvarP, nvarsP, vars)
    return ccall(
        (:GRBgetgenconstrAnd, libgurobi),
        Cint,
        (Ptr{GRBmodel}, Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}),
        model,
        genconstr,
        resvarP,
        nvarsP,
        vars,
    )
end

function GRBgetgenconstrOr(model, genconstr, resvarP, nvarsP, vars)
    return ccall(
        (:GRBgetgenconstrOr, libgurobi),
        Cint,
        (Ptr{GRBmodel}, Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}),
        model,
        genconstr,
        resvarP,
        nvarsP,
        vars,
    )
end

function GRBgetgenconstrNorm(model, genconstr, resvarP, nvarsP, vars, whichP)
    return ccall(
        (:GRBgetgenconstrNorm, libgurobi),
        Cint,
        (Ptr{GRBmodel}, Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}),
        model,
        genconstr,
        resvarP,
        nvarsP,
        vars,
        whichP,
    )
end

function GRBgetgenconstrNL(
    model,
    genconstr,
    resvarP,
    nnodesP,
    opcode,
    data,
    parent,
)
    return ccall(
        (:GRBgetgenconstrNL, libgurobi),
        Cint,
        (
            Ptr{GRBmodel},
            Cint,
            Ptr{Cint},
            Ptr{Cint},
            Ptr{Cint},
            Ptr{Cdouble},
            Ptr{Cint},
        ),
        model,
        genconstr,
        resvarP,
        nnodesP,
        opcode,
        data,
        parent,
    )
end

function GRBgetgenconstrIndicator(
    model,
    genconstr,
    binvarP,
    binvalP,
    nvarsP,
    vars,
    vals,
    senseP,
    rhsP,
)
    return ccall(
        (:GRBgetgenconstrIndicator, libgurobi),
        Cint,
        (
            Ptr{GRBmodel},
            Cint,
            Ptr{Cint},
            Ptr{Cint},
            Ptr{Cint},
            Ptr{Cint},
            Ptr{Cdouble},
            Ptr{Cchar},
            Ptr{Cdouble},
        ),
        model,
        genconstr,
        binvarP,
        binvalP,
        nvarsP,
        vars,
        vals,
        senseP,
        rhsP,
    )
end

function GRBgetgenconstrPWL(model, genconstr, xvarP, yvarP, nptsP, xpts, ypts)
    return ccall(
        (:GRBgetgenconstrPWL, libgurobi),
        Cint,
        (
            Ptr{GRBmodel},
            Cint,
            Ptr{Cint},
            Ptr{Cint},
            Ptr{Cint},
            Ptr{Cdouble},
            Ptr{Cdouble},
        ),
        model,
        genconstr,
        xvarP,
        yvarP,
        nptsP,
        xpts,
        ypts,
    )
end

function GRBgetgenconstrPoly(model, genconstr, xvarP, yvarP, plenP, p)
    return ccall(
        (:GRBgetgenconstrPoly, libgurobi),
        Cint,
        (Ptr{GRBmodel}, Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}),
        model,
        genconstr,
        xvarP,
        yvarP,
        plenP,
        p,
    )
end

function GRBgetgenconstrExpA(model, genconstr, xvarP, yvarP, aP)
    return ccall(
        (:GRBgetgenconstrExpA, libgurobi),
        Cint,
        (Ptr{GRBmodel}, Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}),
        model,
        genconstr,
        xvarP,
        yvarP,
        aP,
    )
end

function GRBgetgenconstrLogA(model, genconstr, xvarP, yvarP, aP)
    return ccall(
        (:GRBgetgenconstrLogA, libgurobi),
        Cint,
        (Ptr{GRBmodel}, Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}),
        model,
        genconstr,
        xvarP,
        yvarP,
        aP,
    )
end

function GRBgetgenconstrPow(model, genconstr, xvarP, yvarP, aP)
    return ccall(
        (:GRBgetgenconstrPow, libgurobi),
        Cint,
        (Ptr{GRBmodel}, Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}),
        model,
        genconstr,
        xvarP,
        yvarP,
        aP,
    )
end

function GRBgetgenconstrExp(model, genconstr, xvarP, yvarP)
    return ccall(
        (:GRBgetgenconstrExp, libgurobi),
        Cint,
        (Ptr{GRBmodel}, Cint, Ptr{Cint}, Ptr{Cint}),
        model,
        genconstr,
        xvarP,
        yvarP,
    )
end

function GRBgetgenconstrLog(model, genconstr, xvarP, yvarP)
    return ccall(
        (:GRBgetgenconstrLog, libgurobi),
        Cint,
        (Ptr{GRBmodel}, Cint, Ptr{Cint}, Ptr{Cint}),
        model,
        genconstr,
        xvarP,
        yvarP,
    )
end

function GRBgetgenconstrLogistic(model, genconstr, xvarP, yvarP)
    return ccall(
        (:GRBgetgenconstrLogistic, libgurobi),
        Cint,
        (Ptr{GRBmodel}, Cint, Ptr{Cint}, Ptr{Cint}),
        model,
        genconstr,
        xvarP,
        yvarP,
    )
end

function GRBgetgenconstrSin(model, genconstr, xvarP, yvarP)
    return ccall(
        (:GRBgetgenconstrSin, libgurobi),
        Cint,
        (Ptr{GRBmodel}, Cint, Ptr{Cint}, Ptr{Cint}),
        model,
        genconstr,
        xvarP,
        yvarP,
    )
end

function GRBgetgenconstrCos(model, genconstr, xvarP, yvarP)
    return ccall(
        (:GRBgetgenconstrCos, libgurobi),
        Cint,
        (Ptr{GRBmodel}, Cint, Ptr{Cint}, Ptr{Cint}),
        model,
        genconstr,
        xvarP,
        yvarP,
    )
end

function GRBgetgenconstrTan(model, genconstr, xvarP, yvarP)
    return ccall(
        (:GRBgetgenconstrTan, libgurobi),
        Cint,
        (Ptr{GRBmodel}, Cint, Ptr{Cint}, Ptr{Cint}),
        model,
        genconstr,
        xvarP,
        yvarP,
    )
end

function GRBgetq(model, numqnzP, qrow, qcol, qval)
    return ccall(
        (:GRBgetq, libgurobi),
        Cint,
        (Ptr{GRBmodel}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}),
        model,
        numqnzP,
        qrow,
        qcol,
        qval,
    )
end

function GRBgetqconstr(
    model,
    qconstr,
    numlnzP,
    lind,
    lval,
    numqnzP,
    qrow,
    qcol,
    qval,
)
    return ccall(
        (:GRBgetqconstr, libgurobi),
        Cint,
        (
            Ptr{GRBmodel},
            Cint,
            Ptr{Cint},
            Ptr{Cint},
            Ptr{Cdouble},
            Ptr{Cint},
            Ptr{Cint},
            Ptr{Cint},
            Ptr{Cdouble},
        ),
        model,
        qconstr,
        numlnzP,
        lind,
        lval,
        numqnzP,
        qrow,
        qcol,
        qval,
    )
end

function GRBgetvarbyname(model, name, indexP)
    return ccall(
        (:GRBgetvarbyname, libgurobi),
        Cint,
        (Ptr{GRBmodel}, Ptr{Cchar}, Ptr{Cint}),
        model,
        name,
        indexP,
    )
end

function GRBgetconstrbyname(model, name, indexP)
    return ccall(
        (:GRBgetconstrbyname, libgurobi),
        Cint,
        (Ptr{GRBmodel}, Ptr{Cchar}, Ptr{Cint}),
        model,
        name,
        indexP,
    )
end

function GRBgetqconstrbyname(model, name, indexP)
    return ccall(
        (:GRBgetqconstrbyname, libgurobi),
        Cint,
        (Ptr{GRBmodel}, Ptr{Cchar}, Ptr{Cint}),
        model,
        name,
        indexP,
    )
end

function GRBgetpwlobj(model, var, pointsP, x, y)
    return ccall(
        (:GRBgetpwlobj, libgurobi),
        Cint,
        (Ptr{GRBmodel}, Cint, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cdouble}),
        model,
        var,
        pointsP,
        x,
        y,
    )
end

function GRBoptimize(model)
    return ccall((:GRBoptimize, libgurobi), Cint, (Ptr{GRBmodel},), model)
end

function GRBoptimizeasync(model)
    return ccall((:GRBoptimizeasync, libgurobi), Cint, (Ptr{GRBmodel},), model)
end

function GRBoptimizebatch(model, batchid)
    return ccall(
        (:GRBoptimizebatch, libgurobi),
        Cint,
        (Ptr{GRBmodel}, Ptr{Cchar}),
        model,
        batchid,
    )
end

function GRBcopymodel(model)
    return ccall(
        (:GRBcopymodel, libgurobi),
        Ptr{GRBmodel},
        (Ptr{GRBmodel},),
        model,
    )
end

function GRBcopymodeltoenv(model, env, resultP)
    return ccall(
        (:GRBcopymodeltoenv, libgurobi),
        Cint,
        (Ptr{GRBmodel}, Ptr{GRBenv}, Ptr{Ptr{GRBmodel}}),
        model,
        env,
        resultP,
    )
end

function GRBfixmodel(model, fixedP)
    return ccall(
        (:GRBfixmodel, libgurobi),
        Cint,
        (Ptr{GRBmodel}, Ptr{Ptr{GRBmodel}}),
        model,
        fixedP,
    )
end

function GRBfeasrelax(
    model,
    relaxobjtype,
    minrelax,
    lbpen,
    ubpen,
    rhspen,
    feasobjP,
)
    return ccall(
        (:GRBfeasrelax, libgurobi),
        Cint,
        (
            Ptr{GRBmodel},
            Cint,
            Cint,
            Ptr{Cdouble},
            Ptr{Cdouble},
            Ptr{Cdouble},
            Ptr{Cdouble},
        ),
        model,
        relaxobjtype,
        minrelax,
        lbpen,
        ubpen,
        rhspen,
        feasobjP,
    )
end

function GRBsinglescenariomodel(model, singlescenarioP)
    return ccall(
        (:GRBsinglescenariomodel, libgurobi),
        Cint,
        (Ptr{GRBmodel}, Ptr{Ptr{GRBmodel}}),
        model,
        singlescenarioP,
    )
end

function GRBconverttofixed(model)
    return ccall((:GRBconverttofixed, libgurobi), Cint, (Ptr{GRBmodel},), model)
end

function GRBgetcbwhatinfo(cbdata, what, typeP, sizeP)
    return ccall(
        (:GRBgetcbwhatinfo, libgurobi),
        Cint,
        (Ptr{Cvoid}, Cint, Ptr{Cint}, Ptr{Cint}),
        cbdata,
        what,
        typeP,
        sizeP,
    )
end

function GRBrelaxmodel(model, relaxedP)
    return ccall(
        (:GRBrelaxmodel, libgurobi),
        Cint,
        (Ptr{GRBmodel}, Ptr{Ptr{GRBmodel}}),
        model,
        relaxedP,
    )
end

function GRBpresolvemodel(model, presolvedP)
    return ccall(
        (:GRBpresolvemodel, libgurobi),
        Cint,
        (Ptr{GRBmodel}, Ptr{Ptr{GRBmodel}}),
        model,
        presolvedP,
    )
end

function GRBiismodel(model, iisP)
    return ccall(
        (:GRBiismodel, libgurobi),
        Cint,
        (Ptr{GRBmodel}, Ptr{Ptr{GRBmodel}}),
        model,
        iisP,
    )
end

function GRBfeasibility(model, feasP)
    return ccall(
        (:GRBfeasibility, libgurobi),
        Cint,
        (Ptr{GRBmodel}, Ptr{Ptr{GRBmodel}}),
        model,
        feasP,
    )
end

function GRBlinearizemodel(model, linearizedP)
    return ccall(
        (:GRBlinearizemodel, libgurobi),
        Cint,
        (Ptr{GRBmodel}, Ptr{Ptr{GRBmodel}}),
        model,
        linearizedP,
    )
end

function GRBresultmodel(model, type, resultP)
    return ccall(
        (:GRBresultmodel, libgurobi),
        Cint,
        (Ptr{GRBmodel}, Ptr{Cchar}, Ptr{Ptr{GRBmodel}}),
        model,
        type,
        resultP,
    )
end

function GRBfixedmodel(model)
    return ccall(
        (:GRBfixedmodel, libgurobi),
        Ptr{GRBmodel},
        (Ptr{GRBmodel},),
        model,
    )
end

function GRBdualmodel(model, dualP)
    return ccall(
        (:GRBdualmodel, libgurobi),
        Cint,
        (Ptr{GRBmodel}, Ptr{Ptr{GRBmodel}}),
        model,
        dualP,
    )
end

function GRBreadmodel(env, filename, modelP)
    return ccall(
        (:GRBreadmodel, libgurobi),
        Cint,
        (Ptr{GRBenv}, Ptr{Cchar}, Ptr{Ptr{GRBmodel}}),
        env,
        filename,
        modelP,
    )
end

function GRBread(model, filename)
    return ccall(
        (:GRBread, libgurobi),
        Cint,
        (Ptr{GRBmodel}, Ptr{Cchar}),
        model,
        filename,
    )
end

function GRBwrite(model, filename)
    return ccall(
        (:GRBwrite, libgurobi),
        Cint,
        (Ptr{GRBmodel}, Ptr{Cchar}),
        model,
        filename,
    )
end

function GRBismodelfile(filename)
    return ccall((:GRBismodelfile, libgurobi), Cint, (Ptr{Cchar},), filename)
end

function GRBisattrfile(filename)
    return ccall((:GRBisattrfile, libgurobi), Cint, (Ptr{Cchar},), filename)
end

function GRBfiletype(filename)
    return ccall((:GRBfiletype, libgurobi), Cint, (Ptr{Cchar},), filename)
end

function GRBisrecordfile(filename)
    return ccall((:GRBisrecordfile, libgurobi), Cint, (Ptr{Cchar},), filename)
end

function GRBgetjsonsolution(model, buffP)
    return ccall(
        (:GRBgetjsonsolution, libgurobi),
        Cint,
        (Ptr{GRBmodel}, Ptr{Ptr{Cchar}}),
        model,
        buffP,
    )
end

function GRBloadjson(env, fname, buffP)
    return ccall(
        (:GRBloadjson, libgurobi),
        Cint,
        (Ptr{GRBenv}, Ptr{Cchar}, Ptr{Ptr{Cchar}}),
        env,
        fname,
        buffP,
    )
end

function GRBnewmodel(env, modelP, Pname, numvars, obj, lb, ub, vtype, varnames)
    return ccall(
        (:GRBnewmodel, libgurobi),
        Cint,
        (
            Ptr{GRBenv},
            Ptr{Ptr{GRBmodel}},
            Ptr{Cchar},
            Cint,
            Ptr{Cdouble},
            Ptr{Cdouble},
            Ptr{Cdouble},
            Ptr{Cchar},
            Ptr{Ptr{Cchar}},
        ),
        env,
        modelP,
        Pname,
        numvars,
        obj,
        lb,
        ub,
        vtype,
        varnames,
    )
end

function GRBloadmodel(
    env,
    modelP,
    Pname,
    numvars,
    numconstrs,
    objsense,
    objcon,
    obj,
    sense,
    rhs,
    vbeg,
    vlen,
    vind,
    vval,
    lb,
    ub,
    vtype,
    varnames,
    constrnames,
)
    return ccall(
        (:GRBloadmodel, libgurobi),
        Cint,
        (
            Ptr{GRBenv},
            Ptr{Ptr{GRBmodel}},
            Ptr{Cchar},
            Cint,
            Cint,
            Cint,
            Cdouble,
            Ptr{Cdouble},
            Ptr{Cchar},
            Ptr{Cdouble},
            Ptr{Cint},
            Ptr{Cint},
            Ptr{Cint},
            Ptr{Cdouble},
            Ptr{Cdouble},
            Ptr{Cdouble},
            Ptr{Cchar},
            Ptr{Ptr{Cchar}},
            Ptr{Ptr{Cchar}},
        ),
        env,
        modelP,
        Pname,
        numvars,
        numconstrs,
        objsense,
        objcon,
        obj,
        sense,
        rhs,
        vbeg,
        vlen,
        vind,
        vval,
        lb,
        ub,
        vtype,
        varnames,
        constrnames,
    )
end

function GRBXloadmodel(
    env,
    modelP,
    Pname,
    numvars,
    numconstrs,
    objsense,
    objcon,
    obj,
    sense,
    rhs,
    vbeg,
    vlen,
    vind,
    vval,
    lb,
    ub,
    vtype,
    varnames,
    constrnames,
)
    return ccall(
        (:GRBXloadmodel, libgurobi),
        Cint,
        (
            Ptr{GRBenv},
            Ptr{Ptr{GRBmodel}},
            Ptr{Cchar},
            Cint,
            Cint,
            Cint,
            Cdouble,
            Ptr{Cdouble},
            Ptr{Cchar},
            Ptr{Cdouble},
            Ptr{Csize_t},
            Ptr{Cint},
            Ptr{Cint},
            Ptr{Cdouble},
            Ptr{Cdouble},
            Ptr{Cdouble},
            Ptr{Cchar},
            Ptr{Ptr{Cchar}},
            Ptr{Ptr{Cchar}},
        ),
        env,
        modelP,
        Pname,
        numvars,
        numconstrs,
        objsense,
        objcon,
        obj,
        sense,
        rhs,
        vbeg,
        vlen,
        vind,
        vval,
        lb,
        ub,
        vtype,
        varnames,
        constrnames,
    )
end

function GRBaddvar(model, numnz, vind, vval, obj, lb, ub, vtype, varname)
    return ccall(
        (:GRBaddvar, libgurobi),
        Cint,
        (
            Ptr{GRBmodel},
            Cint,
            Ptr{Cint},
            Ptr{Cdouble},
            Cdouble,
            Cdouble,
            Cdouble,
            Cchar,
            Ptr{Cchar},
        ),
        model,
        numnz,
        vind,
        vval,
        obj,
        lb,
        ub,
        vtype,
        varname,
    )
end

function GRBaddvars(
    model,
    numvars,
    numnz,
    vbeg,
    vind,
    vval,
    obj,
    lb,
    ub,
    vtype,
    varnames,
)
    return ccall(
        (:GRBaddvars, libgurobi),
        Cint,
        (
            Ptr{GRBmodel},
            Cint,
            Cint,
            Ptr{Cint},
            Ptr{Cint},
            Ptr{Cdouble},
            Ptr{Cdouble},
            Ptr{Cdouble},
            Ptr{Cdouble},
            Ptr{Cchar},
            Ptr{Ptr{Cchar}},
        ),
        model,
        numvars,
        numnz,
        vbeg,
        vind,
        vval,
        obj,
        lb,
        ub,
        vtype,
        varnames,
    )
end

function GRBXaddvars(
    model,
    numvars,
    numnz,
    vbeg,
    vind,
    vval,
    obj,
    lb,
    ub,
    vtype,
    varnames,
)
    return ccall(
        (:GRBXaddvars, libgurobi),
        Cint,
        (
            Ptr{GRBmodel},
            Cint,
            Csize_t,
            Ptr{Csize_t},
            Ptr{Cint},
            Ptr{Cdouble},
            Ptr{Cdouble},
            Ptr{Cdouble},
            Ptr{Cdouble},
            Ptr{Cchar},
            Ptr{Ptr{Cchar}},
        ),
        model,
        numvars,
        numnz,
        vbeg,
        vind,
        vval,
        obj,
        lb,
        ub,
        vtype,
        varnames,
    )
end

function GRBaddconstr(model, numnz, cind, cval, sense, rhs, constrname)
    return ccall(
        (:GRBaddconstr, libgurobi),
        Cint,
        (
            Ptr{GRBmodel},
            Cint,
            Ptr{Cint},
            Ptr{Cdouble},
            Cchar,
            Cdouble,
            Ptr{Cchar},
        ),
        model,
        numnz,
        cind,
        cval,
        sense,
        rhs,
        constrname,
    )
end

function GRBaddconstrs(
    model,
    numconstrs,
    numnz,
    cbeg,
    cind,
    cval,
    sense,
    rhs,
    constrnames,
)
    return ccall(
        (:GRBaddconstrs, libgurobi),
        Cint,
        (
            Ptr{GRBmodel},
            Cint,
            Cint,
            Ptr{Cint},
            Ptr{Cint},
            Ptr{Cdouble},
            Ptr{Cchar},
            Ptr{Cdouble},
            Ptr{Ptr{Cchar}},
        ),
        model,
        numconstrs,
        numnz,
        cbeg,
        cind,
        cval,
        sense,
        rhs,
        constrnames,
    )
end

function GRBXaddconstrs(
    model,
    numconstrs,
    numnz,
    cbeg,
    cind,
    cval,
    sense,
    rhs,
    constrnames,
)
    return ccall(
        (:GRBXaddconstrs, libgurobi),
        Cint,
        (
            Ptr{GRBmodel},
            Cint,
            Csize_t,
            Ptr{Csize_t},
            Ptr{Cint},
            Ptr{Cdouble},
            Ptr{Cchar},
            Ptr{Cdouble},
            Ptr{Ptr{Cchar}},
        ),
        model,
        numconstrs,
        numnz,
        cbeg,
        cind,
        cval,
        sense,
        rhs,
        constrnames,
    )
end

function GRBaddrangeconstr(model, numnz, cind, cval, lower, upper, constrname)
    return ccall(
        (:GRBaddrangeconstr, libgurobi),
        Cint,
        (
            Ptr{GRBmodel},
            Cint,
            Ptr{Cint},
            Ptr{Cdouble},
            Cdouble,
            Cdouble,
            Ptr{Cchar},
        ),
        model,
        numnz,
        cind,
        cval,
        lower,
        upper,
        constrname,
    )
end

function GRBaddrangeconstrs(
    model,
    numconstrs,
    numnz,
    cbeg,
    cind,
    cval,
    lower,
    upper,
    constrnames,
)
    return ccall(
        (:GRBaddrangeconstrs, libgurobi),
        Cint,
        (
            Ptr{GRBmodel},
            Cint,
            Cint,
            Ptr{Cint},
            Ptr{Cint},
            Ptr{Cdouble},
            Ptr{Cdouble},
            Ptr{Cdouble},
            Ptr{Ptr{Cchar}},
        ),
        model,
        numconstrs,
        numnz,
        cbeg,
        cind,
        cval,
        lower,
        upper,
        constrnames,
    )
end

function GRBXaddrangeconstrs(
    model,
    numconstrs,
    numnz,
    cbeg,
    cind,
    cval,
    lower,
    upper,
    constrnames,
)
    return ccall(
        (:GRBXaddrangeconstrs, libgurobi),
        Cint,
        (
            Ptr{GRBmodel},
            Cint,
            Csize_t,
            Ptr{Csize_t},
            Ptr{Cint},
            Ptr{Cdouble},
            Ptr{Cdouble},
            Ptr{Cdouble},
            Ptr{Ptr{Cchar}},
        ),
        model,
        numconstrs,
        numnz,
        cbeg,
        cind,
        cval,
        lower,
        upper,
        constrnames,
    )
end

function GRBaddsos(model, numsos, nummembers, types, beg, ind, weight)
    return ccall(
        (:GRBaddsos, libgurobi),
        Cint,
        (
            Ptr{GRBmodel},
            Cint,
            Cint,
            Ptr{Cint},
            Ptr{Cint},
            Ptr{Cint},
            Ptr{Cdouble},
        ),
        model,
        numsos,
        nummembers,
        types,
        beg,
        ind,
        weight,
    )
end

function GRBaddgenconstrMax(model, name, resvar, nvars, vars, constant)
    return ccall(
        (:GRBaddgenconstrMax, libgurobi),
        Cint,
        (Ptr{GRBmodel}, Ptr{Cchar}, Cint, Cint, Ptr{Cint}, Cdouble),
        model,
        name,
        resvar,
        nvars,
        vars,
        constant,
    )
end

function GRBaddgenconstrMin(model, name, resvar, nvars, vars, constant)
    return ccall(
        (:GRBaddgenconstrMin, libgurobi),
        Cint,
        (Ptr{GRBmodel}, Ptr{Cchar}, Cint, Cint, Ptr{Cint}, Cdouble),
        model,
        name,
        resvar,
        nvars,
        vars,
        constant,
    )
end

function GRBaddgenconstrAbs(model, name, resvar, argvar)
    return ccall(
        (:GRBaddgenconstrAbs, libgurobi),
        Cint,
        (Ptr{GRBmodel}, Ptr{Cchar}, Cint, Cint),
        model,
        name,
        resvar,
        argvar,
    )
end

function GRBaddgenconstrAnd(model, name, resvar, nvars, vars)
    return ccall(
        (:GRBaddgenconstrAnd, libgurobi),
        Cint,
        (Ptr{GRBmodel}, Ptr{Cchar}, Cint, Cint, Ptr{Cint}),
        model,
        name,
        resvar,
        nvars,
        vars,
    )
end

function GRBaddgenconstrOr(model, name, resvar, nvars, vars)
    return ccall(
        (:GRBaddgenconstrOr, libgurobi),
        Cint,
        (Ptr{GRBmodel}, Ptr{Cchar}, Cint, Cint, Ptr{Cint}),
        model,
        name,
        resvar,
        nvars,
        vars,
    )
end

function GRBaddgenconstrNorm(model, name, resvar, nvars, vars, which)
    return ccall(
        (:GRBaddgenconstrNorm, libgurobi),
        Cint,
        (Ptr{GRBmodel}, Ptr{Cchar}, Cint, Cint, Ptr{Cint}, Cdouble),
        model,
        name,
        resvar,
        nvars,
        vars,
        which,
    )
end

function GRBaddgenconstrNL(model, name, resvar, nnodes, opcode, data, parent)
    return ccall(
        (:GRBaddgenconstrNL, libgurobi),
        Cint,
        (
            Ptr{GRBmodel},
            Ptr{Cchar},
            Cint,
            Cint,
            Ptr{Cint},
            Ptr{Cdouble},
            Ptr{Cint},
        ),
        model,
        name,
        resvar,
        nnodes,
        opcode,
        data,
        parent,
    )
end

function GRBaddgenconstrIndicator(
    model,
    name,
    binvar,
    binval,
    nvars,
    vars,
    vals,
    sense,
    rhs,
)
    return ccall(
        (:GRBaddgenconstrIndicator, libgurobi),
        Cint,
        (
            Ptr{GRBmodel},
            Ptr{Cchar},
            Cint,
            Cint,
            Cint,
            Ptr{Cint},
            Ptr{Cdouble},
            Cchar,
            Cdouble,
        ),
        model,
        name,
        binvar,
        binval,
        nvars,
        vars,
        vals,
        sense,
        rhs,
    )
end

function GRBaddgenconstrPWL(model, name, xvar, yvar, npts, xpts, ypts)
    return ccall(
        (:GRBaddgenconstrPWL, libgurobi),
        Cint,
        (
            Ptr{GRBmodel},
            Ptr{Cchar},
            Cint,
            Cint,
            Cint,
            Ptr{Cdouble},
            Ptr{Cdouble},
        ),
        model,
        name,
        xvar,
        yvar,
        npts,
        xpts,
        ypts,
    )
end

function GRBaddgenconstrPoly(model, name, xvar, yvar, plen, p, options)
    return ccall(
        (:GRBaddgenconstrPoly, libgurobi),
        Cint,
        (Ptr{GRBmodel}, Ptr{Cchar}, Cint, Cint, Cint, Ptr{Cdouble}, Ptr{Cchar}),
        model,
        name,
        xvar,
        yvar,
        plen,
        p,
        options,
    )
end

function GRBaddgenconstrExpA(model, name, xvar, yvar, a, options)
    return ccall(
        (:GRBaddgenconstrExpA, libgurobi),
        Cint,
        (Ptr{GRBmodel}, Ptr{Cchar}, Cint, Cint, Cdouble, Ptr{Cchar}),
        model,
        name,
        xvar,
        yvar,
        a,
        options,
    )
end

function GRBaddgenconstrLogA(model, name, xvar, yvar, a, options)
    return ccall(
        (:GRBaddgenconstrLogA, libgurobi),
        Cint,
        (Ptr{GRBmodel}, Ptr{Cchar}, Cint, Cint, Cdouble, Ptr{Cchar}),
        model,
        name,
        xvar,
        yvar,
        a,
        options,
    )
end

function GRBaddgenconstrPow(model, name, xvar, yvar, a, options)
    return ccall(
        (:GRBaddgenconstrPow, libgurobi),
        Cint,
        (Ptr{GRBmodel}, Ptr{Cchar}, Cint, Cint, Cdouble, Ptr{Cchar}),
        model,
        name,
        xvar,
        yvar,
        a,
        options,
    )
end

function GRBaddgenconstrExp(model, name, xvar, yvar, options)
    return ccall(
        (:GRBaddgenconstrExp, libgurobi),
        Cint,
        (Ptr{GRBmodel}, Ptr{Cchar}, Cint, Cint, Ptr{Cchar}),
        model,
        name,
        xvar,
        yvar,
        options,
    )
end

function GRBaddgenconstrLog(model, name, xvar, yvar, options)
    return ccall(
        (:GRBaddgenconstrLog, libgurobi),
        Cint,
        (Ptr{GRBmodel}, Ptr{Cchar}, Cint, Cint, Ptr{Cchar}),
        model,
        name,
        xvar,
        yvar,
        options,
    )
end

function GRBaddgenconstrSin(model, name, xvar, yvar, options)
    return ccall(
        (:GRBaddgenconstrSin, libgurobi),
        Cint,
        (Ptr{GRBmodel}, Ptr{Cchar}, Cint, Cint, Ptr{Cchar}),
        model,
        name,
        xvar,
        yvar,
        options,
    )
end

function GRBaddgenconstrCos(model, name, xvar, yvar, options)
    return ccall(
        (:GRBaddgenconstrCos, libgurobi),
        Cint,
        (Ptr{GRBmodel}, Ptr{Cchar}, Cint, Cint, Ptr{Cchar}),
        model,
        name,
        xvar,
        yvar,
        options,
    )
end

function GRBaddgenconstrTan(model, name, xvar, yvar, options)
    return ccall(
        (:GRBaddgenconstrTan, libgurobi),
        Cint,
        (Ptr{GRBmodel}, Ptr{Cchar}, Cint, Cint, Ptr{Cchar}),
        model,
        name,
        xvar,
        yvar,
        options,
    )
end

function GRBaddgenconstrLogistic(model, name, xvar, yvar, options)
    return ccall(
        (:GRBaddgenconstrLogistic, libgurobi),
        Cint,
        (Ptr{GRBmodel}, Ptr{Cchar}, Cint, Cint, Ptr{Cchar}),
        model,
        name,
        xvar,
        yvar,
        options,
    )
end

function GRBaddqconstr(
    model,
    numlnz,
    lind,
    lval,
    numqnz,
    qrow,
    qcol,
    qval,
    sense,
    rhs,
    QCname,
)
    return ccall(
        (:GRBaddqconstr, libgurobi),
        Cint,
        (
            Ptr{GRBmodel},
            Cint,
            Ptr{Cint},
            Ptr{Cdouble},
            Cint,
            Ptr{Cint},
            Ptr{Cint},
            Ptr{Cdouble},
            Cchar,
            Cdouble,
            Ptr{Cchar},
        ),
        model,
        numlnz,
        lind,
        lval,
        numqnz,
        qrow,
        qcol,
        qval,
        sense,
        rhs,
        QCname,
    )
end

function GRBaddcone(model, nummembers, members)
    return ccall(
        (:GRBaddcone, libgurobi),
        Cint,
        (Ptr{GRBmodel}, Cint, Ptr{Cint}),
        model,
        nummembers,
        members,
    )
end

function GRBaddqpterms(model, numqnz, qrow, qcol, qval)
    return ccall(
        (:GRBaddqpterms, libgurobi),
        Cint,
        (Ptr{GRBmodel}, Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}),
        model,
        numqnz,
        qrow,
        qcol,
        qval,
    )
end

function GRBdelvars(model, len, ind)
    return ccall(
        (:GRBdelvars, libgurobi),
        Cint,
        (Ptr{GRBmodel}, Cint, Ptr{Cint}),
        model,
        len,
        ind,
    )
end

function GRBdelconstrs(model, len, ind)
    return ccall(
        (:GRBdelconstrs, libgurobi),
        Cint,
        (Ptr{GRBmodel}, Cint, Ptr{Cint}),
        model,
        len,
        ind,
    )
end

function GRBdelsos(model, len, ind)
    return ccall(
        (:GRBdelsos, libgurobi),
        Cint,
        (Ptr{GRBmodel}, Cint, Ptr{Cint}),
        model,
        len,
        ind,
    )
end

function GRBdelgenconstrs(model, len, ind)
    return ccall(
        (:GRBdelgenconstrs, libgurobi),
        Cint,
        (Ptr{GRBmodel}, Cint, Ptr{Cint}),
        model,
        len,
        ind,
    )
end

function GRBdelqconstrs(model, len, ind)
    return ccall(
        (:GRBdelqconstrs, libgurobi),
        Cint,
        (Ptr{GRBmodel}, Cint, Ptr{Cint}),
        model,
        len,
        ind,
    )
end

function GRBdelq(model)
    return ccall((:GRBdelq, libgurobi), Cint, (Ptr{GRBmodel},), model)
end

function GRBchgcoeffs(model, cnt, cind, vind, val)
    return ccall(
        (:GRBchgcoeffs, libgurobi),
        Cint,
        (Ptr{GRBmodel}, Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}),
        model,
        cnt,
        cind,
        vind,
        val,
    )
end

function GRBXchgcoeffs(model, cnt, cind, vind, val)
    return ccall(
        (:GRBXchgcoeffs, libgurobi),
        Cint,
        (Ptr{GRBmodel}, Csize_t, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}),
        model,
        cnt,
        cind,
        vind,
        val,
    )
end

function GRBsetpwlobj(model, var, points, x, y)
    return ccall(
        (:GRBsetpwlobj, libgurobi),
        Cint,
        (Ptr{GRBmodel}, Cint, Cint, Ptr{Cdouble}, Ptr{Cdouble}),
        model,
        var,
        points,
        x,
        y,
    )
end

function GRBupdatemodel(model)
    return ccall((:GRBupdatemodel, libgurobi), Cint, (Ptr{GRBmodel},), model)
end

function GRBreset(model, clearall)
    return ccall(
        (:GRBreset, libgurobi),
        Cint,
        (Ptr{GRBmodel}, Cint),
        model,
        clearall,
    )
end

function GRBfreemodel(model)
    return ccall((:GRBfreemodel, libgurobi), Cint, (Ptr{GRBmodel},), model)
end

function GRBcomputeIIS(model)
    return ccall((:GRBcomputeIIS, libgurobi), Cint, (Ptr{GRBmodel},), model)
end

struct _GRBsvec
    len::Cint
    ind::Ptr{Cint}
    val::Ptr{Cdouble}
end

const GRBsvec = _GRBsvec

function GRBFSolve(model, b, x)
    return ccall(
        (:GRBFSolve, libgurobi),
        Cint,
        (Ptr{GRBmodel}, Ptr{GRBsvec}, Ptr{GRBsvec}),
        model,
        b,
        x,
    )
end

function GRBBinvColj(model, j, x)
    return ccall(
        (:GRBBinvColj, libgurobi),
        Cint,
        (Ptr{GRBmodel}, Cint, Ptr{GRBsvec}),
        model,
        j,
        x,
    )
end

function GRBBinvj(model, j, x)
    return ccall(
        (:GRBBinvj, libgurobi),
        Cint,
        (Ptr{GRBmodel}, Cint, Ptr{GRBsvec}),
        model,
        j,
        x,
    )
end

function GRBBSolve(model, b, x)
    return ccall(
        (:GRBBSolve, libgurobi),
        Cint,
        (Ptr{GRBmodel}, Ptr{GRBsvec}, Ptr{GRBsvec}),
        model,
        b,
        x,
    )
end

function GRBBinvi(model, i, x)
    return ccall(
        (:GRBBinvi, libgurobi),
        Cint,
        (Ptr{GRBmodel}, Cint, Ptr{GRBsvec}),
        model,
        i,
        x,
    )
end

function GRBBinvRowi(model, i, x)
    return ccall(
        (:GRBBinvRowi, libgurobi),
        Cint,
        (Ptr{GRBmodel}, Cint, Ptr{GRBsvec}),
        model,
        i,
        x,
    )
end

function GRBgetBasisHead(model, bhead)
    return ccall(
        (:GRBgetBasisHead, libgurobi),
        Cint,
        (Ptr{GRBmodel}, Ptr{Cint}),
        model,
        bhead,
    )
end

function GRBcbstoponemultiobj(model, cbdata, objnum)
    return ccall(
        (:GRBcbstoponemultiobj, libgurobi),
        Cint,
        (Ptr{GRBmodel}, Ptr{Cvoid}, Cint),
        model,
        cbdata,
        objnum,
    )
end

function GRBsingularvectors(model, left, right)
    return ccall(
        (:GRBsingularvectors, libgurobi),
        Cint,
        (Ptr{GRBmodel}, Ptr{Cdouble}, Ptr{Cdouble}),
        model,
        left,
        right,
    )
end

function GRBstrongbranch(model, num, cand, downobjbd, upobjbd, statusP)
    return ccall(
        (:GRBstrongbranch, libgurobi),
        Cint,
        (Ptr{GRBmodel}, Cint, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint}),
        model,
        num,
        cand,
        downobjbd,
        upobjbd,
        statusP,
    )
end

function GRBcheckmodel(model)
    return ccall((:GRBcheckmodel, libgurobi), Cint, (Ptr{GRBmodel},), model)
end

function GRBterminate(model)
    return ccall((:GRBterminate, libgurobi), Cvoid, (Ptr{GRBmodel},), model)
end

function GRBreplay(filename)
    return ccall((:GRBreplay, libgurobi), Cint, (Ptr{Cchar},), filename)
end

function GRBsetobjective(
    model,
    sense,
    constant,
    lnz,
    lind,
    lval,
    qnz,
    qrow,
    qcol,
    qval,
)
    return ccall(
        (:GRBsetobjective, libgurobi),
        Cint,
        (
            Ptr{GRBmodel},
            Cint,
            Cdouble,
            Cint,
            Ptr{Cint},
            Ptr{Cdouble},
            Cint,
            Ptr{Cint},
            Ptr{Cint},
            Ptr{Cdouble},
        ),
        model,
        sense,
        constant,
        lnz,
        lind,
        lval,
        qnz,
        qrow,
        qcol,
        qval,
    )
end

function GRBsetobjectiven(
    model,
    index,
    priority,
    weight,
    abstol,
    reltol,
    name,
    constant,
    lnz,
    lind,
    lval,
)
    return ccall(
        (:GRBsetobjectiven, libgurobi),
        Cint,
        (
            Ptr{GRBmodel},
            Cint,
            Cint,
            Cdouble,
            Cdouble,
            Cdouble,
            Ptr{Cchar},
            Cdouble,
            Cint,
            Ptr{Cint},
            Ptr{Cdouble},
        ),
        model,
        index,
        priority,
        weight,
        abstol,
        reltol,
        name,
        constant,
        lnz,
        lind,
        lval,
    )
end

function GRBclean2(lenP, ind, val)
    return ccall(
        (:GRBclean2, libgurobi),
        Cvoid,
        (Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}),
        lenP,
        ind,
        val,
    )
end

function GRBclean3(lenP, ind0, ind1, val)
    return ccall(
        (:GRBclean3, libgurobi),
        Cvoid,
        (Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}),
        lenP,
        ind0,
        ind1,
        val,
    )
end

function GRBprintquality(model)
    return ccall((:GRBprintquality, libgurobi), Cint, (Ptr{GRBmodel},), model)
end

function GRBmsg(env, message)
    return ccall(
        (:GRBmsg, libgurobi),
        Cvoid,
        (Ptr{GRBenv}, Ptr{Cchar}),
        env,
        message,
    )
end

function GRBgetintparam(env, paramname, valueP)
    return ccall(
        (:GRBgetintparam, libgurobi),
        Cint,
        (Ptr{GRBenv}, Ptr{Cchar}, Ptr{Cint}),
        env,
        paramname,
        valueP,
    )
end

function GRBgetdblparam(env, paramname, valueP)
    return ccall(
        (:GRBgetdblparam, libgurobi),
        Cint,
        (Ptr{GRBenv}, Ptr{Cchar}, Ptr{Cdouble}),
        env,
        paramname,
        valueP,
    )
end

function GRBgetstrparam(env, paramname, valueP)
    return ccall(
        (:GRBgetstrparam, libgurobi),
        Cint,
        (Ptr{GRBenv}, Ptr{Cchar}, Ptr{Cchar}),
        env,
        paramname,
        valueP,
    )
end

function GRBgetlongstrparam(env, paramname, valueP, size, requiredlenP)
    return ccall(
        (:GRBgetlongstrparam, libgurobi),
        Cint,
        (Ptr{GRBenv}, Ptr{Cchar}, Ptr{Cchar}, Cint, Ptr{Cint}),
        env,
        paramname,
        valueP,
        size,
        requiredlenP,
    )
end

function GRBgetintparaminfo(env, paramname, valueP, minP, maxP, defP)
    return ccall(
        (:GRBgetintparaminfo, libgurobi),
        Cint,
        (Ptr{GRBenv}, Ptr{Cchar}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}),
        env,
        paramname,
        valueP,
        minP,
        maxP,
        defP,
    )
end

function GRBgetdblparaminfo(env, paramname, valueP, minP, maxP, defP)
    return ccall(
        (:GRBgetdblparaminfo, libgurobi),
        Cint,
        (
            Ptr{GRBenv},
            Ptr{Cchar},
            Ptr{Cdouble},
            Ptr{Cdouble},
            Ptr{Cdouble},
            Ptr{Cdouble},
        ),
        env,
        paramname,
        valueP,
        minP,
        maxP,
        defP,
    )
end

function GRBgetstrparaminfo(env, paramname, valueP, defP)
    return ccall(
        (:GRBgetstrparaminfo, libgurobi),
        Cint,
        (Ptr{GRBenv}, Ptr{Cchar}, Ptr{Cchar}, Ptr{Cchar}),
        env,
        paramname,
        valueP,
        defP,
    )
end

function GRBgetparamflags(env, parname, valueP)
    return ccall(
        (:GRBgetparamflags, libgurobi),
        Cint,
        (Ptr{GRBenv}, Ptr{Cchar}, Ptr{Cuint}),
        env,
        parname,
        valueP,
    )
end

function GRBsetparam(env, paramname, value)
    return ccall(
        (:GRBsetparam, libgurobi),
        Cint,
        (Ptr{GRBenv}, Ptr{Cchar}, Ptr{Cchar}),
        env,
        paramname,
        value,
    )
end

function GRBsetintparam(env, paramname, value)
    return ccall(
        (:GRBsetintparam, libgurobi),
        Cint,
        (Ptr{GRBenv}, Ptr{Cchar}, Cint),
        env,
        paramname,
        value,
    )
end

function GRBsetdblparam(env, paramname, value)
    return ccall(
        (:GRBsetdblparam, libgurobi),
        Cint,
        (Ptr{GRBenv}, Ptr{Cchar}, Cdouble),
        env,
        paramname,
        value,
    )
end

function GRBsetstrparam(env, paramname, value)
    return ccall(
        (:GRBsetstrparam, libgurobi),
        Cint,
        (Ptr{GRBenv}, Ptr{Cchar}, Ptr{Cchar}),
        env,
        paramname,
        value,
    )
end

function GRBgetparamtype(env, paramname)
    return ccall(
        (:GRBgetparamtype, libgurobi),
        Cint,
        (Ptr{GRBenv}, Ptr{Cchar}),
        env,
        paramname,
    )
end

function GRBresetparams(env)
    return ccall((:GRBresetparams, libgurobi), Cint, (Ptr{GRBenv},), env)
end

function GRBcopyparams(dest, src)
    return ccall(
        (:GRBcopyparams, libgurobi),
        Cint,
        (Ptr{GRBenv}, Ptr{GRBenv}),
        dest,
        src,
    )
end

function GRBwriteparams(env, filename)
    return ccall(
        (:GRBwriteparams, libgurobi),
        Cint,
        (Ptr{GRBenv}, Ptr{Cchar}),
        env,
        filename,
    )
end

function GRBreadparams(env, filename)
    return ccall(
        (:GRBreadparams, libgurobi),
        Cint,
        (Ptr{GRBenv}, Ptr{Cchar}),
        env,
        filename,
    )
end

function GRBreadconcurrentsettings(model, filename)
    return ccall(
        (:GRBreadconcurrentsettings, libgurobi),
        Cint,
        (Ptr{GRBmodel}, Ptr{Cchar}),
        model,
        filename,
    )
end

function GRBreadmultiobjsettings(model, filename)
    return ccall(
        (:GRBreadmultiobjsettings, libgurobi),
        Cint,
        (Ptr{GRBmodel}, Ptr{Cchar}),
        model,
        filename,
    )
end

function GRBreadtunebasesettings(env, filename)
    return ccall(
        (:GRBreadtunebasesettings, libgurobi),
        Cint,
        (Ptr{GRBenv}, Ptr{Cchar}),
        env,
        filename,
    )
end

function GRBgetnumparams(env)
    return ccall((:GRBgetnumparams, libgurobi), Cint, (Ptr{GRBenv},), env)
end

function GRBgetparamname(env, parnum, paramnameP)
    return ccall(
        (:GRBgetparamname, libgurobi),
        Cint,
        (Ptr{GRBenv}, Cint, Ptr{Ptr{Cchar}}),
        env,
        parnum,
        paramnameP,
    )
end

function GRBgetnumattributes(model)
    return ccall(
        (:GRBgetnumattributes, libgurobi),
        Cint,
        (Ptr{GRBmodel},),
        model,
    )
end

function GRBgetattrname(model, i, attrnameP)
    return ccall(
        (:GRBgetattrname, libgurobi),
        Cint,
        (Ptr{GRBmodel}, Cint, Ptr{Ptr{Cchar}}),
        model,
        i,
        attrnameP,
    )
end

function GRBstartenv(env)
    return ccall((:GRBstartenv, libgurobi), Cint, (Ptr{GRBenv},), env)
end

function GRBloadenvadv(
    envP,
    logfilename,
    apitype,
    major,
    minor,
    tech,
    server,
    router,
    password,
    group,
    priority,
    idletimeout,
    cloudaccessid,
    cloudsecretkey,
    cb,
    usrdata,
    logcb,
    logdata,
)
    return ccall(
        (:GRBloadenvadv, libgurobi),
        Cint,
        (
            Ptr{Ptr{GRBenv}},
            Ptr{Cchar},
            Cint,
            Cint,
            Cint,
            Cint,
            Ptr{Cchar},
            Ptr{Cchar},
            Ptr{Cchar},
            Ptr{Cchar},
            Cint,
            Cint,
            Ptr{Cchar},
            Ptr{Cchar},
            Ptr{Cvoid},
            Ptr{Cvoid},
            Ptr{Cvoid},
            Ptr{Cvoid},
        ),
        envP,
        logfilename,
        apitype,
        major,
        minor,
        tech,
        server,
        router,
        password,
        group,
        priority,
        idletimeout,
        cloudaccessid,
        cloudsecretkey,
        cb,
        usrdata,
        logcb,
        logdata,
    )
end

function GRBgetenv(model)
    return ccall((:GRBgetenv, libgurobi), Ptr{GRBenv}, (Ptr{GRBmodel},), model)
end

function GRBgetconcurrentenv(model, num)
    return ccall(
        (:GRBgetconcurrentenv, libgurobi),
        Ptr{GRBenv},
        (Ptr{GRBmodel}, Cint),
        model,
        num,
    )
end

function GRBdiscardconcurrentenvs(model)
    return ccall(
        (:GRBdiscardconcurrentenvs, libgurobi),
        Cvoid,
        (Ptr{GRBmodel},),
        model,
    )
end

function GRBgetmultiobjenv(model, num)
    return ccall(
        (:GRBgetmultiobjenv, libgurobi),
        Ptr{GRBenv},
        (Ptr{GRBmodel}, Cint),
        model,
        num,
    )
end

function GRBdiscardmultiobjenvs(model)
    return ccall(
        (:GRBdiscardmultiobjenvs, libgurobi),
        Cvoid,
        (Ptr{GRBmodel},),
        model,
    )
end

function GRBgettuneenv(env, num)
    return ccall(
        (:GRBgettuneenv, libgurobi),
        Ptr{GRBenv},
        (Ptr{GRBenv}, Cint),
        env,
        num,
    )
end

function GRBdiscardtuneenvs(env)
    return ccall((:GRBdiscardtuneenvs, libgurobi), Cvoid, (Ptr{GRBenv},), env)
end

function GRBreleaselicense(env)
    return ccall((:GRBreleaselicense, libgurobi), Cvoid, (Ptr{GRBenv},), env)
end

function GRBfreeenv(env)
    return ccall((:GRBfreeenv, libgurobi), Cvoid, (Ptr{GRBenv},), env)
end

function GRBgeterrormsg(env)
    return ccall((:GRBgeterrormsg, libgurobi), Ptr{Cchar}, (Ptr{GRBenv},), env)
end

function GRBgetmerrormsg(model)
    return ccall(
        (:GRBgetmerrormsg, libgurobi),
        Ptr{Cchar},
        (Ptr{GRBmodel},),
        model,
    )
end

function GRBgetcommstats(
    env,
    recvtimeP,
    recvbytesP,
    recvmsgsP,
    sendtimeP,
    sendbytesP,
    sendmsgsP,
)
    return ccall(
        (:GRBgetcommstats, libgurobi),
        Cvoid,
        (
            Ptr{GRBenv},
            Ptr{Cdouble},
            Ptr{Cdouble},
            Ptr{Cdouble},
            Ptr{Cdouble},
            Ptr{Cdouble},
            Ptr{Cdouble},
        ),
        env,
        recvtimeP,
        recvbytesP,
        recvmsgsP,
        sendtimeP,
        sendbytesP,
        sendmsgsP,
    )
end

function GRBversion(majorP, minorP, technicalP)
    return ccall(
        (:GRBversion, libgurobi),
        Cvoid,
        (Ptr{Cint}, Ptr{Cint}, Ptr{Cint}),
        majorP,
        minorP,
        technicalP,
    )
end

function GRBgetdistro(str)
    return ccall((:GRBgetdistro, libgurobi), Cvoid, (Ptr{Cchar},), str)
end

function GRBplatform()
    return ccall((:GRBplatform, libgurobi), Ptr{Cchar}, ())
end

function GRBplatformext()
    return ccall((:GRBplatformext, libgurobi), Ptr{Cchar}, ())
end

function GRBlisttokens()
    return ccall((:GRBlisttokens, libgurobi), Cint, ())
end

function GRBgetwlstokenlifespan(env, lifespanP)
    return ccall(
        (:GRBgetwlstokenlifespan, libgurobi),
        Cint,
        (Ptr{GRBenv}, Ptr{Cint}),
        env,
        lifespanP,
    )
end

function GRBsortIDi(len, ind, val)
    return ccall(
        (:GRBsortIDi, libgurobi),
        Cvoid,
        (Cint, Ptr{Cint}, Ptr{Cdouble}),
        len,
        ind,
        val,
    )
end

function GRBabortbatch(batch)
    return ccall((:GRBabortbatch, libgurobi), Cint, (Ptr{GRBbatch},), batch)
end

function GRBdiscardbatch(batch)
    return ccall((:GRBdiscardbatch, libgurobi), Cint, (Ptr{GRBbatch},), batch)
end

function GRBretrybatch(batch)
    return ccall((:GRBretrybatch, libgurobi), Cint, (Ptr{GRBbatch},), batch)
end

function GRBfreebatch(batch)
    return ccall((:GRBfreebatch, libgurobi), Cint, (Ptr{GRBbatch},), batch)
end

function GRBgetbatch(env, batchID, batchP)
    return ccall(
        (:GRBgetbatch, libgurobi),
        Cint,
        (Ptr{GRBenv}, Ptr{Cchar}, Ptr{Ptr{GRBbatch}}),
        env,
        batchID,
        batchP,
    )
end

function GRBgetbatchjsonsolution(batch, jsonsolP)
    return ccall(
        (:GRBgetbatchjsonsolution, libgurobi),
        Cint,
        (Ptr{GRBbatch}, Ptr{Ptr{Cchar}}),
        batch,
        jsonsolP,
    )
end

function GRBgetbatchintattr(batch, attrname, valueP)
    return ccall(
        (:GRBgetbatchintattr, libgurobi),
        Cint,
        (Ptr{GRBbatch}, Ptr{Cchar}, Ptr{Cint}),
        batch,
        attrname,
        valueP,
    )
end

function GRBgetbatchstrattr(batch, attrname, valueP)
    return ccall(
        (:GRBgetbatchstrattr, libgurobi),
        Cint,
        (Ptr{GRBbatch}, Ptr{Cchar}, Ptr{Ptr{Cchar}}),
        batch,
        attrname,
        valueP,
    )
end

function GRBgetbatchattrname(env, n, attrnameP)
    return ccall(
        (:GRBgetbatchattrname, libgurobi),
        Cint,
        (Ptr{GRBenv}, Cint, Ptr{Ptr{Cchar}}),
        env,
        n,
        attrnameP,
    )
end

function GRBgetbatchattrflags(batch, attrname, flagsP)
    return ccall(
        (:GRBgetbatchattrflags, libgurobi),
        Cint,
        (Ptr{GRBbatch}, Ptr{Cchar}, Ptr{Cuint}),
        batch,
        attrname,
        flagsP,
    )
end

function GRBgetbatchattrinfo(batch, attrname, datatypeP, settableP)
    return ccall(
        (:GRBgetbatchattrinfo, libgurobi),
        Cint,
        (Ptr{GRBbatch}, Ptr{Cchar}, Ptr{Cint}, Ptr{Cint}),
        batch,
        attrname,
        datatypeP,
        settableP,
    )
end

function GRBupdatebatch(batch)
    return ccall((:GRBupdatebatch, libgurobi), Cint, (Ptr{GRBbatch},), batch)
end

function GRBwritebatchjsonsolution(batch, filename)
    return ccall(
        (:GRBwritebatchjsonsolution, libgurobi),
        Cint,
        (Ptr{GRBbatch}, Ptr{Cchar}),
        batch,
        filename,
    )
end

function GRBgetnumbatchattributes(env)
    return ccall(
        (:GRBgetnumbatchattributes, libgurobi),
        Cint,
        (Ptr{GRBenv},),
        env,
    )
end

function GRBgetbatchenv(batch)
    return ccall(
        (:GRBgetbatchenv, libgurobi),
        Ptr{GRBenv},
        (Ptr{GRBbatch},),
        batch,
    )
end

function GRBfree(ptr)
    return ccall((:GRBfree, libgurobi), Cvoid, (Ptr{Cvoid},), ptr)
end

function GRBsync(model)
    return ccall((:GRBsync, libgurobi), Cint, (Ptr{GRBmodel},), model)
end

function GRBpingserver(server, password)
    return ccall(
        (:GRBpingserver, libgurobi),
        Cint,
        (Ptr{Cchar}, Ptr{Cchar}),
        server,
        password,
    )
end

function GRBprefetchattr(model, attrname)
    return ccall(
        (:GRBprefetchattr, libgurobi),
        Cint,
        (Ptr{GRBmodel}, Ptr{Cchar}),
        model,
        attrname,
    )
end

function GRBtunemodel(model)
    return ccall((:GRBtunemodel, libgurobi), Cint, (Ptr{GRBmodel},), model)
end

function GRBtunemodels(env, nummodels, models)
    return ccall(
        (:GRBtunemodels, libgurobi),
        Cint,
        (Ptr{GRBenv}, Cint, Ptr{Ptr{GRBmodel}}),
        env,
        nummodels,
        models,
    )
end

function GRBgettuneresult(model, i)
    return ccall(
        (:GRBgettuneresult, libgurobi),
        Cint,
        (Ptr{GRBmodel}, Cint),
        model,
        i,
    )
end

function GRBgettunelog(model, i, logP)
    return ccall(
        (:GRBgettunelog, libgurobi),
        Cint,
        (Ptr{GRBmodel}, Cint, Ptr{Ptr{Cchar}}),
        model,
        i,
        logP,
    )
end

function GRBwritetunelog(model, result, filename)
    return ccall(
        (:GRBwritetunelog, libgurobi),
        Cint,
        (Ptr{GRBmodel}, Cint, Ptr{Cchar}),
        model,
        result,
        filename,
    )
end

function GRBtuneparamsPrint()
    return ccall((:GRBtuneparamsPrint, libgurobi), Cvoid, ())
end

# Skipping MacroDefinition: __cdecl

# Skipping MacroDefinition: GRB_DEPRECATED ( WHY , X ) X __attribute__ ( ( deprecated ) )

const GRB_VERSION_MAJOR = 12

const GRB_VERSION_MINOR = 0

const GRB_VERSION_TECHNICAL = 0

const DEFAULT_CS_PRIORITY = 0

const MAX_CS_PRIORITY = 100

const DEFAULT_CS_PORT = 61000

const DEFAULT_CS_HANGUP = 60

const GRB_C_MIN_ERROR = 10001

const GRB_ERROR_OUT_OF_MEMORY = 10001

const GRB_ERROR_NULL_ARGUMENT = 10002

const GRB_ERROR_INVALID_ARGUMENT = 10003

const GRB_ERROR_UNKNOWN_ATTRIBUTE = 10004

const GRB_ERROR_DATA_NOT_AVAILABLE = 10005

const GRB_ERROR_INDEX_OUT_OF_RANGE = 10006

const GRB_ERROR_UNKNOWN_PARAMETER = 10007

const GRB_ERROR_VALUE_OUT_OF_RANGE = 10008

const GRB_ERROR_NO_LICENSE = 10009

const GRB_ERROR_SIZE_LIMIT_EXCEEDED = 10010

const GRB_ERROR_CALLBACK = 10011

const GRB_ERROR_FILE_READ = 10012

const GRB_ERROR_FILE_WRITE = 10013

const GRB_ERROR_NUMERIC = 10014

const GRB_ERROR_IIS_NOT_INFEASIBLE = 10015

const GRB_ERROR_NOT_FOR_MIP = 10016

const GRB_ERROR_OPTIMIZATION_IN_PROGRESS = 10017

const GRB_ERROR_DUPLICATES = 10018

const GRB_ERROR_NODEFILE = 10019

const GRB_ERROR_Q_NOT_PSD = 10020

const GRB_ERROR_QCP_EQUALITY_CONSTRAINT = 10021

const GRB_ERROR_NETWORK = 10022

const GRB_ERROR_JOB_REJECTED = 10023

const GRB_ERROR_NOT_SUPPORTED = 10024

const GRB_ERROR_EXCEED_2B_NONZEROS = 10025

const GRB_ERROR_INVALID_PIECEWISE_OBJ = 10026

const GRB_ERROR_UPDATEMODE_CHANGE = 10027

const GRB_ERROR_CLOUD = 10028

const GRB_ERROR_MODEL_MODIFICATION = 10029

const GRB_ERROR_CSWORKER = 10030

const GRB_ERROR_TUNE_MODEL_TYPES = 10031

const GRB_ERROR_SECURITY = 10032

const GRB_C_MAX_ERROR = 10032

const GRB_MINIMIZE = 1

const GRB_MAXIMIZE = -1

const GRB_SOS_TYPE1 = 1

const GRB_SOS_TYPE2 = 2

const GRB_INFINITY = 1.0e100

const GRB_UNDEFINED = 1.0e101

const GRB_MAXINT = 2000000000

const GRB_MAX_NAMELEN = 255

const GRB_MAX_STRLEN = 512

const GRB_MAX_TAGLEN = 10240

const GRB_MAX_CONCURRENT = 64

const GRB_INT_ATTR_NUMCONSTRS = "NumConstrs"

const GRB_INT_ATTR_NUMVARS = "NumVars"

const GRB_INT_ATTR_NUMSOS = "NumSOS"

const GRB_INT_ATTR_NUMQCONSTRS = "NumQConstrs"

const GRB_INT_ATTR_NUMGENCONSTRS = "NumGenConstrs"

const GRB_INT_ATTR_NUMNZS = "NumNZs"

const GRB_DBL_ATTR_DNUMNZS = "DNumNZs"

const GRB_INT_ATTR_NUMQNZS = "NumQNZs"

const GRB_INT_ATTR_NUMQCNZS = "NumQCNZs"

const GRB_INT_ATTR_NUMINTVARS = "NumIntVars"

const GRB_INT_ATTR_NUMBINVARS = "NumBinVars"

const GRB_INT_ATTR_NUMPWLOBJVARS = "NumPWLObjVars"

const GRB_STR_ATTR_MODELNAME = "ModelName"

const GRB_INT_ATTR_MODELSENSE = "ModelSense"

const GRB_DBL_ATTR_OBJCON = "ObjCon"

const GRB_INT_ATTR_IS_MIP = "IsMIP"

const GRB_INT_ATTR_IS_QP = "IsQP"

const GRB_INT_ATTR_IS_QCP = "IsQCP"

const GRB_INT_ATTR_IS_MULTIOBJ = "IsMultiObj"

const GRB_INT_ATTR_LICENSE_EXPIRATION = "LicenseExpiration"

const GRB_INT_ATTR_NUMTAGGED = "NumTagged"

const GRB_INT_ATTR_FINGERPRINT = "Fingerprint"

const GRB_INT_ATTR_BATCHERRORCODE = "BatchErrorCode"

const GRB_STR_ATTR_BATCHERRORMESSAGE = "BatchErrorMessage"

const GRB_STR_ATTR_BATCHID = "BatchID"

const GRB_INT_ATTR_BATCHSTATUS = "BatchStatus"

const GRB_DBL_ATTR_LB = "LB"

const GRB_DBL_ATTR_UB = "UB"

const GRB_DBL_ATTR_OBJ = "Obj"

const GRB_CHAR_ATTR_VTYPE = "VType"

const GRB_DBL_ATTR_START = "Start"

const GRB_DBL_ATTR_PSTART = "PStart"

const GRB_INT_ATTR_BRANCHPRIORITY = "BranchPriority"

const GRB_STR_ATTR_VARNAME = "VarName"

const GRB_INT_ATTR_PWLOBJCVX = "PWLObjCvx"

const GRB_DBL_ATTR_VARHINTVAL = "VarHintVal"

const GRB_INT_ATTR_VARHINTPRI = "VarHintPri"

const GRB_INT_ATTR_PARTITION = "Partition"

const GRB_INT_ATTR_POOLIGNORE = "PoolIgnore"

const GRB_STR_ATTR_VTAG = "VTag"

const GRB_STR_ATTR_CTAG = "CTag"

const GRB_DBL_ATTR_RHS = "RHS"

const GRB_DBL_ATTR_DSTART = "DStart"

const GRB_CHAR_ATTR_SENSE = "Sense"

const GRB_STR_ATTR_CONSTRNAME = "ConstrName"

const GRB_INT_ATTR_LAZY = "Lazy"

const GRB_STR_ATTR_QCTAG = "QCTag"

const GRB_DBL_ATTR_QCRHS = "QCRHS"

const GRB_CHAR_ATTR_QCSENSE = "QCSense"

const GRB_STR_ATTR_QCNAME = "QCName"

const GRB_INT_ATTR_GENCONSTRTYPE = "GenConstrType"

const GRB_STR_ATTR_GENCONSTRNAME = "GenConstrName"

const GRB_INT_ATTR_FUNCPIECES = "FuncPieces"

const GRB_DBL_ATTR_FUNCPIECEERROR = "FuncPieceError"

const GRB_DBL_ATTR_FUNCPIECELENGTH = "FuncPieceLength"

const GRB_DBL_ATTR_FUNCPIECERATIO = "FuncPieceRatio"

const GRB_INT_ATTR_FUNCNONLINEAR = "FuncNonlinear"

const GRB_DBL_ATTR_MAX_COEFF = "MaxCoeff"

const GRB_DBL_ATTR_MIN_COEFF = "MinCoeff"

const GRB_DBL_ATTR_MAX_BOUND = "MaxBound"

const GRB_DBL_ATTR_MIN_BOUND = "MinBound"

const GRB_DBL_ATTR_MAX_OBJ_COEFF = "MaxObjCoeff"

const GRB_DBL_ATTR_MIN_OBJ_COEFF = "MinObjCoeff"

const GRB_DBL_ATTR_MAX_RHS = "MaxRHS"

const GRB_DBL_ATTR_MIN_RHS = "MinRHS"

const GRB_DBL_ATTR_MAX_QCCOEFF = "MaxQCCoeff"

const GRB_DBL_ATTR_MIN_QCCOEFF = "MinQCCoeff"

const GRB_DBL_ATTR_MAX_QOBJ_COEFF = "MaxQObjCoeff"

const GRB_DBL_ATTR_MIN_QOBJ_COEFF = "MinQObjCoeff"

const GRB_DBL_ATTR_MAX_QCLCOEFF = "MaxQCLCoeff"

const GRB_DBL_ATTR_MIN_QCLCOEFF = "MinQCLCoeff"

const GRB_DBL_ATTR_MAX_QCRHS = "MaxQCRHS"

const GRB_DBL_ATTR_MIN_QCRHS = "MinQCRHS"

const GRB_DBL_ATTR_RUNTIME = "Runtime"

const GRB_DBL_ATTR_WORK = "Work"

const GRB_INT_ATTR_STATUS = "Status"

const GRB_DBL_ATTR_OBJVAL = "ObjVal"

const GRB_DBL_ATTR_OBJBOUND = "ObjBound"

const GRB_DBL_ATTR_OBJBOUNDC = "ObjBoundC"

const GRB_DBL_ATTR_POOLOBJBOUND = "PoolObjBound"

const GRB_DBL_ATTR_POOLOBJVAL = "PoolObjVal"

const GRB_DBL_ATTR_MIPGAP = "MIPGap"

const GRB_INT_ATTR_SOLCOUNT = "SolCount"

const GRB_DBL_ATTR_ITERCOUNT = "IterCount"

const GRB_INT_ATTR_BARITERCOUNT = "BarIterCount"

const GRB_DBL_ATTR_NODECOUNT = "NodeCount"

const GRB_DBL_ATTR_OPENNODECOUNT = "OpenNodeCount"

const GRB_INT_ATTR_HASDUALNORM = "HasDualNorm"

const GRB_INT_ATTR_CONCURRENTWINMETHOD = "ConcurrentWinMethod"

const GRB_DBL_ATTR_X = "X"

const GRB_DBL_ATTR_XN = "Xn"

const GRB_DBL_ATTR_BARX = "BarX"

const GRB_DBL_ATTR_BARPI = "BarPi"

const GRB_DBL_ATTR_RC = "RC"

const GRB_DBL_ATTR_VDUALNORM = "VDualNorm"

const GRB_INT_ATTR_VBASIS = "VBasis"

const GRB_DBL_ATTR_PI = "Pi"

const GRB_DBL_ATTR_QCPI = "QCPi"

const GRB_DBL_ATTR_SLACK = "Slack"

const GRB_DBL_ATTR_QCSLACK = "QCSlack"

const GRB_DBL_ATTR_CDUALNORM = "CDualNorm"

const GRB_INT_ATTR_CBASIS = "CBasis"

const GRB_DBL_ATTR_MAX_VIO = "MaxVio"

const GRB_DBL_ATTR_BOUND_VIO = "BoundVio"

const GRB_DBL_ATTR_BOUND_SVIO = "BoundSVio"

const GRB_INT_ATTR_BOUND_VIO_INDEX = "BoundVioIndex"

const GRB_INT_ATTR_BOUND_SVIO_INDEX = "BoundSVioIndex"

const GRB_DBL_ATTR_BOUND_VIO_SUM = "BoundVioSum"

const GRB_DBL_ATTR_BOUND_SVIO_SUM = "BoundSVioSum"

const GRB_DBL_ATTR_CONSTR_VIO = "ConstrVio"

const GRB_DBL_ATTR_CONSTR_SVIO = "ConstrSVio"

const GRB_INT_ATTR_CONSTR_VIO_INDEX = "ConstrVioIndex"

const GRB_INT_ATTR_CONSTR_SVIO_INDEX = "ConstrSVioIndex"

const GRB_DBL_ATTR_CONSTR_VIO_SUM = "ConstrVioSum"

const GRB_DBL_ATTR_CONSTR_SVIO_SUM = "ConstrSVioSum"

const GRB_DBL_ATTR_CONSTR_RESIDUAL = "ConstrResidual"

const GRB_DBL_ATTR_CONSTR_SRESIDUAL = "ConstrSResidual"

const GRB_INT_ATTR_CONSTR_RESIDUAL_INDEX = "ConstrResidualIndex"

const GRB_INT_ATTR_CONSTR_SRESIDUAL_INDEX = "ConstrSResidualIndex"

const GRB_DBL_ATTR_CONSTR_RESIDUAL_SUM = "ConstrResidualSum"

const GRB_DBL_ATTR_CONSTR_SRESIDUAL_SUM = "ConstrSResidualSum"

const GRB_DBL_ATTR_DUAL_VIO = "DualVio"

const GRB_DBL_ATTR_DUAL_SVIO = "DualSVio"

const GRB_INT_ATTR_DUAL_VIO_INDEX = "DualVioIndex"

const GRB_INT_ATTR_DUAL_SVIO_INDEX = "DualSVioIndex"

const GRB_DBL_ATTR_DUAL_VIO_SUM = "DualVioSum"

const GRB_DBL_ATTR_DUAL_SVIO_SUM = "DualSVioSum"

const GRB_DBL_ATTR_DUAL_RESIDUAL = "DualResidual"

const GRB_DBL_ATTR_DUAL_SRESIDUAL = "DualSResidual"

const GRB_INT_ATTR_DUAL_RESIDUAL_INDEX = "DualResidualIndex"

const GRB_INT_ATTR_DUAL_SRESIDUAL_INDEX = "DualSResidualIndex"

const GRB_DBL_ATTR_DUAL_RESIDUAL_SUM = "DualResidualSum"

const GRB_DBL_ATTR_DUAL_SRESIDUAL_SUM = "DualSResidualSum"

const GRB_DBL_ATTR_INT_VIO = "IntVio"

const GRB_INT_ATTR_INT_VIO_INDEX = "IntVioIndex"

const GRB_DBL_ATTR_INT_VIO_SUM = "IntVioSum"

const GRB_DBL_ATTR_COMPL_VIO = "ComplVio"

const GRB_INT_ATTR_COMPL_VIO_INDEX = "ComplVioIndex"

const GRB_DBL_ATTR_COMPL_VIO_SUM = "ComplVioSum"

const GRB_DBL_ATTR_KAPPA = "Kappa"

const GRB_DBL_ATTR_KAPPA_EXACT = "KappaExact"

const GRB_DBL_ATTR_N2KAPPA = "N2Kappa"

const GRB_DBL_ATTR_SA_OBJLOW = "SAObjLow"

const GRB_DBL_ATTR_SA_OBJUP = "SAObjUp"

const GRB_DBL_ATTR_SA_LBLOW = "SALBLow"

const GRB_DBL_ATTR_SA_LBUP = "SALBUp"

const GRB_DBL_ATTR_SA_UBLOW = "SAUBLow"

const GRB_DBL_ATTR_SA_UBUP = "SAUBUp"

const GRB_DBL_ATTR_SA_RHSLOW = "SARHSLow"

const GRB_DBL_ATTR_SA_RHSUP = "SARHSUp"

const GRB_INT_ATTR_IIS_MINIMAL = "IISMinimal"

const GRB_INT_ATTR_IIS_LB = "IISLB"

const GRB_INT_ATTR_IIS_UB = "IISUB"

const GRB_INT_ATTR_IIS_CONSTR = "IISConstr"

const GRB_INT_ATTR_IIS_SOS = "IISSOS"

const GRB_INT_ATTR_IIS_QCONSTR = "IISQConstr"

const GRB_INT_ATTR_IIS_GENCONSTR = "IISGenConstr"

const GRB_INT_ATTR_IIS_LBFORCE = "IISLBForce"

const GRB_INT_ATTR_IIS_UBFORCE = "IISUBForce"

const GRB_INT_ATTR_IIS_CONSTRFORCE = "IISConstrForce"

const GRB_INT_ATTR_IIS_SOSFORCE = "IISSOSForce"

const GRB_INT_ATTR_IIS_QCONSTRFORCE = "IISQConstrForce"

const GRB_INT_ATTR_IIS_GENCONSTRFORCE = "IISGenConstrForce"

const GRB_INT_ATTR_TUNE_RESULTCOUNT = "TuneResultCount"

const GRB_DBL_ATTR_FARKASDUAL = "FarkasDual"

const GRB_DBL_ATTR_FARKASPROOF = "FarkasProof"

const GRB_DBL_ATTR_UNBDRAY = "UnbdRay"

const GRB_INT_ATTR_INFEASVAR = "InfeasVar"

const GRB_INT_ATTR_UNBDVAR = "UnbdVar"

const GRB_INT_ATTR_VARPRESTAT = "VarPreStat"

const GRB_DBL_ATTR_PREFIXVAL = "PreFixVal"

const GRB_DBL_ATTR_OBJN = "ObjN"

const GRB_DBL_ATTR_OBJNVAL = "ObjNVal"

const GRB_DBL_ATTR_OBJNCON = "ObjNCon"

const GRB_DBL_ATTR_OBJNWEIGHT = "ObjNWeight"

const GRB_INT_ATTR_OBJNPRIORITY = "ObjNPriority"

const GRB_DBL_ATTR_OBJNRELTOL = "ObjNRelTol"

const GRB_DBL_ATTR_OBJNABSTOL = "ObjNAbsTol"

const GRB_STR_ATTR_OBJNNAME = "ObjNName"

const GRB_DBL_ATTR_SCENNLB = "ScenNLB"

const GRB_DBL_ATTR_SCENNUB = "ScenNUB"

const GRB_DBL_ATTR_SCENNOBJ = "ScenNObj"

const GRB_DBL_ATTR_SCENNRHS = "ScenNRHS"

const GRB_STR_ATTR_SCENNNAME = "ScenNName"

const GRB_DBL_ATTR_SCENNX = "ScenNX"

const GRB_DBL_ATTR_SCENNOBJBOUND = "ScenNObjBound"

const GRB_DBL_ATTR_SCENNOBJVAL = "ScenNObjVal"

const GRB_INT_ATTR_NUMOBJ = "NumObj"

const GRB_INT_ATTR_NUMSCENARIOS = "NumScenarios"

const GRB_INT_ATTR_NUMSTART = "NumStart"

const GRB_DBL_ATTR_MEMUSED = "MemUsed"

const GRB_DBL_ATTR_MAXMEMUSED = "MaxMemUsed"

const GRB_DBL_ATTR_Xn = "Xn"

const GRB_GENCONSTR_MAX = 0

const GRB_GENCONSTR_MIN = 1

const GRB_GENCONSTR_ABS = 2

const GRB_GENCONSTR_AND = 3

const GRB_GENCONSTR_OR = 4

const GRB_GENCONSTR_NORM = 5

const GRB_GENCONSTR_NL = 6

const GRB_GENCONSTR_INDICATOR = 7

const GRB_GENCONSTR_PWL = 8

const GRB_GENCONSTR_POLY = 9

const GRB_GENCONSTR_EXP = 10

const GRB_GENCONSTR_EXPA = 11

const GRB_GENCONSTR_LOG = 12

const GRB_GENCONSTR_LOGA = 13

const GRB_GENCONSTR_POW = 14

const GRB_GENCONSTR_SIN = 15

const GRB_GENCONSTR_COS = 16

const GRB_GENCONSTR_TAN = 17

const GRB_GENCONSTR_LOGISTIC = 18

const NUMGENCONSTYPES = 19

const GRB_OPCODE_CONSTANT = 0

const GRB_OPCODE_VARIABLE = 1

const GRB_OPCODE_PLUS = 2

const GRB_OPCODE_MINUS = 3

const GRB_OPCODE_MULTIPLY = 4

const GRB_OPCODE_DIVIDE = 5

const GRB_OPCODE_UMINUS = 6

const GRB_OPCODE_SQUARE = 7

const GRB_OPCODE_SQRT = 8

const GRB_OPCODE_SIN = 9

const GRB_OPCODE_COS = 10

const GRB_OPCODE_TAN = 11

const GRB_OPCODE_POW = 12

const GRB_OPCODE_EXP = 13

const GRB_OPCODE_LOG = 14

const GRB_OPCODE_LOG2 = 15

const GRB_OPCODE_LOG10 = 16

const GRB_OPCODE_LOGISTIC = 17

const GRB_CB_POLLING = 0

const GRB_CB_PRESOLVE = 1

const GRB_CB_SIMPLEX = 2

const GRB_CB_MIP = 3

const GRB_CB_MIPSOL = 4

const GRB_CB_MIPNODE = 5

const GRB_CB_MESSAGE = 6

const GRB_CB_BARRIER = 7

const GRB_CB_MULTIOBJ = 8

const GRB_CB_IIS = 9

const GRB_CB_PRE_COLDEL = 1000

const GRB_CB_PRE_ROWDEL = 1001

const GRB_CB_PRE_SENCHG = 1002

const GRB_CB_PRE_BNDCHG = 1003

const GRB_CB_PRE_COECHG = 1004

const GRB_CB_SPX_ITRCNT = 2000

const GRB_CB_SPX_OBJVAL = 2001

const GRB_CB_SPX_PRIMINF = 2002

const GRB_CB_SPX_DUALINF = 2003

const GRB_CB_SPX_ISPERT = 2004

const GRB_CB_MIP_OBJBST = 3000

const GRB_CB_MIP_OBJBND = 3001

const GRB_CB_MIP_NODCNT = 3002

const GRB_CB_MIP_SOLCNT = 3003

const GRB_CB_MIP_CUTCNT = 3004

const GRB_CB_MIP_NODLFT = 3005

const GRB_CB_MIP_ITRCNT = 3006

const GRB_CB_MIP_OPENSCENARIOS = 3007

const GRB_CB_MIP_PHASE = 3008

const GRB_CB_MIPSOL_SOL = 4001

const GRB_CB_MIPSOL_OBJ = 4002

const GRB_CB_MIPSOL_OBJBST = 4003

const GRB_CB_MIPSOL_OBJBND = 4004

const GRB_CB_MIPSOL_NODCNT = 4005

const GRB_CB_MIPSOL_SOLCNT = 4006

const GRB_CB_MIPSOL_OPENSCENARIOS = 4007

const GRB_CB_MIPSOL_PHASE = 4008

const GRB_CB_MIPNODE_STATUS = 5001

const GRB_CB_MIPNODE_REL = 5002

const GRB_CB_MIPNODE_OBJBST = 5003

const GRB_CB_MIPNODE_OBJBND = 5004

const GRB_CB_MIPNODE_NODCNT = 5005

const GRB_CB_MIPNODE_SOLCNT = 5006

const GRB_CB_MIPNODE_BRVAR = 5007

const GRB_CB_MIPNODE_OPENSCENARIOS = 5008

const GRB_CB_MIPNODE_PHASE = 5009

const GRB_CB_MSG_STRING = 6001

const GRB_CB_RUNTIME = 6002

const GRB_CB_WORK = 6003

const GRB_CB_MEMUSED = 6004

const GRB_CB_MAXMEMUSED = 6005

const GRB_CB_BARRIER_ITRCNT = 7001

const GRB_CB_BARRIER_PRIMOBJ = 7002

const GRB_CB_BARRIER_DUALOBJ = 7003

const GRB_CB_BARRIER_PRIMINF = 7004

const GRB_CB_BARRIER_DUALINF = 7005

const GRB_CB_BARRIER_COMPL = 7006

const GRB_CB_MULTIOBJ_OBJCNT = 8001

const GRB_CB_MULTIOBJ_SOLCNT = 8002

const GRB_CB_MULTIOBJ_SOL = 8003

const GRB_CB_MULTIOBJ_ITRCNT = 8004

const GRB_CB_MULTIOBJ_OBJBST = 8005

const GRB_CB_MULTIOBJ_OBJBND = 8006

const GRB_CB_MULTIOBJ_STATUS = 8007

const GRB_CB_MULTIOBJ_MIPGAP = 8008

const GRB_CB_MULTIOBJ_NODCNT = 8009

const GRB_CB_MULTIOBJ_NODLFT = 8010

const GRB_CB_MULTIOBJ_RUNTIME = 8011

const GRB_CB_MULTIOBJ_WORK = 8012

const GRB_CB_IIS_CONSTRMIN = 9001

const GRB_CB_IIS_CONSTRMAX = 9002

const GRB_CB_IIS_CONSTRGUESS = 9003

const GRB_CB_IIS_BOUNDMIN = 9004

const GRB_CB_IIS_BOUNDMAX = 9005

const GRB_CB_IIS_BOUNDGUESS = 9006

const GRB_FEASRELAX_LINEAR = 0

const GRB_FEASRELAX_QUADRATIC = 1

const GRB_FEASRELAX_CARDINALITY = 2

function GRBemptyenvadv(
    envP,
    malloccb,
    callocbc,
    realloccb,
    freecb,
    threadcreatecb,
    threadjoincb,
    syscbusrdata,
)
    return GRBemptyenvadvinternal(
        envP,
        -1,
        GRB_VERSION_MAJOR,
        GRB_VERSION_MINOR,
        GRB_VERSION_TECHNICAL,
        malloccb,
        callocbc,
        realloccb,
        freecb,
        threadcreatecb,
        threadjoincb,
        syscbusrdata,
    )
end

const GRB_LOADED = 1

const GRB_OPTIMAL = 2

const GRB_INFEASIBLE = 3

const GRB_INF_OR_UNBD = 4

const GRB_UNBOUNDED = 5

const GRB_CUTOFF = 6

const GRB_ITERATION_LIMIT = 7

const GRB_NODE_LIMIT = 8

const GRB_TIME_LIMIT = 9

const GRB_SOLUTION_LIMIT = 10

const GRB_INTERRUPTED = 11

const GRB_NUMERIC = 12

const GRB_SUBOPTIMAL = 13

const GRB_INPROGRESS = 14

const GRB_USER_OBJ_LIMIT = 15

const GRB_WORK_LIMIT = 16

const GRB_MEM_LIMIT = 17

const GRB_BASIC = 0

const GRB_NONBASIC_LOWER = -1

const GRB_NONBASIC_UPPER = -2

const GRB_SUPERBASIC = -3

const GRB_INT_PAR_BARITERLIMIT = "BarIterLimit"

const GRB_DBL_PAR_CUTOFF = "Cutoff"

const GRB_DBL_PAR_ITERATIONLIMIT = "IterationLimit"

const GRB_DBL_PAR_NODELIMIT = "NodeLimit"

const GRB_INT_PAR_SOLUTIONLIMIT = "SolutionLimit"

const GRB_DBL_PAR_TIMELIMIT = "TimeLimit"

const GRB_DBL_PAR_WORKLIMIT = "WorkLimit"

const GRB_DBL_PAR_MEMLIMIT = "MemLimit"

const GRB_DBL_PAR_SOFTMEMLIMIT = "SoftMemLimit"

const GRB_DBL_PAR_BESTOBJSTOP = "BestObjStop"

const GRB_DBL_PAR_BESTBDSTOP = "BestBdStop"

const GRB_DBL_PAR_FEASIBILITYTOL = "FeasibilityTol"

const GRB_DBL_PAR_INTFEASTOL = "IntFeasTol"

const GRB_DBL_PAR_MARKOWITZTOL = "MarkowitzTol"

const GRB_DBL_PAR_MIPGAP = "MIPGap"

const GRB_DBL_PAR_MIPGAPABS = "MIPGapAbs"

const GRB_DBL_PAR_OPTIMALITYTOL = "OptimalityTol"

const GRB_DBL_PAR_PSDTOL = "PSDTol"

const GRB_INT_PAR_METHOD = "Method"

const GRB_INT_PAR_CONCURRENTMETHOD = "ConcurrentMethod"

const GRB_DBL_PAR_PERTURBVALUE = "PerturbValue"

const GRB_DBL_PAR_OBJSCALE = "ObjScale"

const GRB_INT_PAR_SCALEFLAG = "ScaleFlag"

const GRB_INT_PAR_SIMPLEXPRICING = "SimplexPricing"

const GRB_INT_PAR_QUAD = "Quad"

const GRB_INT_PAR_NORMADJUST = "NormAdjust"

const GRB_INT_PAR_SIFTING = "Sifting"

const GRB_INT_PAR_SIFTMETHOD = "SiftMethod"

const GRB_INT_PAR_LPWARMSTART = "LPWarmStart"

const GRB_INT_PAR_NETWORKALG = "NetworkAlg"

const GRB_DBL_PAR_BARCONVTOL = "BarConvTol"

const GRB_INT_PAR_BARCORRECTORS = "BarCorrectors"

const GRB_INT_PAR_BARHOMOGENEOUS = "BarHomogeneous"

const GRB_INT_PAR_BARORDER = "BarOrder"

const GRB_DBL_PAR_BARQCPCONVTOL = "BarQCPConvTol"

const GRB_INT_PAR_CROSSOVER = "Crossover"

const GRB_INT_PAR_CROSSOVERBASIS = "CrossoverBasis"

const GRB_INT_PAR_BRANCHDIR = "BranchDir"

const GRB_INT_PAR_DEGENMOVES = "DegenMoves"

const GRB_INT_PAR_DISCONNECTED = "Disconnected"

const GRB_DBL_PAR_HEURISTICS = "Heuristics"

const GRB_DBL_PAR_IMPROVESTARTGAP = "ImproveStartGap"

const GRB_DBL_PAR_IMPROVESTARTTIME = "ImproveStartTime"

const GRB_DBL_PAR_IMPROVESTARTNODES = "ImproveStartNodes"

const GRB_INT_PAR_INTEGRALITYFOCUS = "IntegralityFocus"

const GRB_INT_PAR_MINRELNODES = "MinRelNodes"

const GRB_INT_PAR_MIPFOCUS = "MIPFocus"

const GRB_INT_PAR_NLPHEUR = "NLPHeur"

const GRB_STR_PAR_NODEFILEDIR = "NodefileDir"

const GRB_DBL_PAR_NODEFILESTART = "NodefileStart"

const GRB_INT_PAR_NODEMETHOD = "NodeMethod"

const GRB_DBL_PAR_NORELHEURTIME = "NoRelHeurTime"

const GRB_DBL_PAR_NORELHEURWORK = "NoRelHeurWork"

const GRB_INT_PAR_OBBT = "OBBT"

const GRB_INT_PAR_PUMPPASSES = "PumpPasses"

const GRB_INT_PAR_RINS = "RINS"

const GRB_STR_PAR_SOLFILES = "SolFiles"

const GRB_INT_PAR_STARTNODELIMIT = "StartNodeLimit"

const GRB_INT_PAR_SUBMIPNODES = "SubMIPNodes"

const GRB_INT_PAR_SYMMETRY = "Symmetry"

const GRB_INT_PAR_VARBRANCH = "VarBranch"

const GRB_INT_PAR_SOLUTIONNUMBER = "SolutionNumber"

const GRB_INT_PAR_ZEROOBJNODES = "ZeroObjNodes"

const GRB_INT_PAR_CUTS = "Cuts"

const GRB_INT_PAR_CLIQUECUTS = "CliqueCuts"

const GRB_INT_PAR_COVERCUTS = "CoverCuts"

const GRB_INT_PAR_FLOWCOVERCUTS = "FlowCoverCuts"

const GRB_INT_PAR_FLOWPATHCUTS = "FlowPathCuts"

const GRB_INT_PAR_GUBCOVERCUTS = "GUBCoverCuts"

const GRB_INT_PAR_IMPLIEDCUTS = "ImpliedCuts"

const GRB_INT_PAR_PROJIMPLIEDCUTS = "ProjImpliedCuts"

const GRB_INT_PAR_MIPSEPCUTS = "MIPSepCuts"

const GRB_INT_PAR_MIRCUTS = "MIRCuts"

const GRB_INT_PAR_STRONGCGCUTS = "StrongCGCuts"

const GRB_INT_PAR_MODKCUTS = "ModKCuts"

const GRB_INT_PAR_ZEROHALFCUTS = "ZeroHalfCuts"

const GRB_INT_PAR_NETWORKCUTS = "NetworkCuts"

const GRB_INT_PAR_SUBMIPCUTS = "SubMIPCuts"

const GRB_INT_PAR_INFPROOFCUTS = "InfProofCuts"

const GRB_INT_PAR_RLTCUTS = "RLTCuts"

const GRB_INT_PAR_RELAXLIFTCUTS = "RelaxLiftCuts"

const GRB_INT_PAR_BQPCUTS = "BQPCuts"

const GRB_INT_PAR_PSDCUTS = "PSDCuts"

const GRB_INT_PAR_LIFTPROJECTCUTS = "LiftProjectCuts"

const GRB_INT_PAR_MIXINGCUTS = "MixingCuts"

const GRB_INT_PAR_DUALIMPLIEDCUTS = "DualImpliedCuts"

const GRB_INT_PAR_CUTAGGPASSES = "CutAggPasses"

const GRB_INT_PAR_CUTPASSES = "CutPasses"

const GRB_INT_PAR_GOMORYPASSES = "GomoryPasses"

const GRB_STR_PAR_WORKERPOOL = "WorkerPool"

const GRB_STR_PAR_WORKERPASSWORD = "WorkerPassword"

const GRB_STR_PAR_COMPUTESERVER = "ComputeServer"

const GRB_STR_PAR_TOKENSERVER = "TokenServer"

const GRB_STR_PAR_SERVERPASSWORD = "ServerPassword"

const GRB_INT_PAR_SERVERTIMEOUT = "ServerTimeout"

const GRB_STR_PAR_CSROUTER = "CSRouter"

const GRB_STR_PAR_CSGROUP = "CSGroup"

const GRB_DBL_PAR_CSQUEUETIMEOUT = "CSQueueTimeout"

const GRB_INT_PAR_CSPRIORITY = "CSPriority"

const GRB_INT_PAR_CSIDLETIMEOUT = "CSIdleTimeout"

const GRB_INT_PAR_CSTLSINSECURE = "CSTLSInsecure"

const GRB_INT_PAR_TSPORT = "TSPort"

const GRB_STR_PAR_CLOUDACCESSID = "CloudAccessID"

const GRB_STR_PAR_CLOUDSECRETKEY = "CloudSecretKey"

const GRB_STR_PAR_CLOUDPOOL = "CloudPool"

const GRB_STR_PAR_CLOUDHOST = "CloudHost"

const GRB_STR_PAR_CSMANAGER = "CSManager"

const GRB_STR_PAR_CSAUTHTOKEN = "CSAuthToken"

const GRB_STR_PAR_CSAPIACCESSID = "CSAPIAccessID"

const GRB_STR_PAR_CSAPISECRET = "CSAPISecret"

const GRB_INT_PAR_CSBATCHMODE = "CSBatchMode"

const GRB_STR_PAR_USERNAME = "Username"

const GRB_STR_PAR_CSAPPNAME = "CSAppName"

const GRB_INT_PAR_CSCLIENTLOG = "CSClientLog"

const GRB_STR_PAR_WLSACCESSID = "WLSAccessID"

const GRB_STR_PAR_WLSSECRET = "WLSSecret"

const GRB_INT_PAR_WLSTOKENDURATION = "WLSTokenDuration"

const GRB_DBL_PAR_WLSTOKENREFRESH = "WLSTokenRefresh"

const GRB_STR_PAR_WLSTOKEN = "WLSToken"

const GRB_INT_PAR_LICENSEID = "LicenseID"

const GRB_STR_PAR_WLSPROXY = "WLSProxy"

const GRB_STR_PAR_WLSCONFIG = "WLSConfig"

const GRB_INT_PAR_AGGREGATE = "Aggregate"

const GRB_INT_PAR_AGGFILL = "AggFill"

const GRB_INT_PAR_CONCURRENTMIP = "ConcurrentMIP"

const GRB_INT_PAR_CONCURRENTJOBS = "ConcurrentJobs"

const GRB_INT_PAR_DISPLAYINTERVAL = "DisplayInterval"

const GRB_INT_PAR_DISTRIBUTEDMIPJOBS = "DistributedMIPJobs"

const GRB_INT_PAR_DUALREDUCTIONS = "DualReductions"

const GRB_DBL_PAR_FEASRELAXBIGM = "FeasRelaxBigM"

const GRB_INT_PAR_IISMETHOD = "IISMethod"

const GRB_INT_PAR_INFUNBDINFO = "InfUnbdInfo"

const GRB_INT_PAR_JSONSOLDETAIL = "JSONSolDetail"

const GRB_INT_PAR_LAZYCONSTRAINTS = "LazyConstraints"

const GRB_STR_PAR_LOGFILE = "LogFile"

const GRB_INT_PAR_LOGTOCONSOLE = "LogToConsole"

const GRB_INT_PAR_MIQCPMETHOD = "MIQCPMethod"

const GRB_INT_PAR_NONCONVEX = "NonConvex"

const GRB_INT_PAR_NUMERICFOCUS = "NumericFocus"

const GRB_INT_PAR_OUTPUTFLAG = "OutputFlag"

const GRB_INT_PAR_PRECRUSH = "PreCrush"

const GRB_INT_PAR_PREDEPROW = "PreDepRow"

const GRB_INT_PAR_PREDUAL = "PreDual"

const GRB_INT_PAR_PREPASSES = "PrePasses"

const GRB_INT_PAR_PREQLINEARIZE = "PreQLinearize"

const GRB_INT_PAR_PRESOLVE = "Presolve"

const GRB_DBL_PAR_PRESOS1BIGM = "PreSOS1BigM"

const GRB_DBL_PAR_PRESOS2BIGM = "PreSOS2BigM"

const GRB_INT_PAR_PRESOS1ENCODING = "PreSOS1Encoding"

const GRB_INT_PAR_PRESOS2ENCODING = "PreSOS2Encoding"

const GRB_INT_PAR_PRESPARSIFY = "PreSparsify"

const GRB_INT_PAR_PREMIQCPFORM = "PreMIQCPForm"

const GRB_INT_PAR_QCPDUAL = "QCPDual"

const GRB_INT_PAR_RECORD = "Record"

const GRB_STR_PAR_RESULTFILE = "ResultFile"

const GRB_INT_PAR_SEED = "Seed"

const GRB_INT_PAR_SOLUTIONTARGET = "SolutionTarget"

const GRB_INT_PAR_THREADS = "Threads"

const GRB_INT_PAR_THREADLIMIT = "ThreadLimit"

const GRB_DBL_PAR_TUNETIMELIMIT = "TuneTimeLimit"

const GRB_INT_PAR_TUNERESULTS = "TuneResults"

const GRB_INT_PAR_TUNECRITERION = "TuneCriterion"

const GRB_INT_PAR_TUNETRIALS = "TuneTrials"

const GRB_INT_PAR_TUNEOUTPUT = "TuneOutput"

const GRB_INT_PAR_TUNEJOBS = "TuneJobs"

const GRB_DBL_PAR_TUNECLEANUP = "TuneCleanup"

const GRB_DBL_PAR_TUNETARGETMIPGAP = "TuneTargetMIPGap"

const GRB_DBL_PAR_TUNETARGETTIME = "TuneTargetTime"

const GRB_INT_PAR_TUNEMETRIC = "TuneMetric"

const GRB_INT_PAR_TUNEDYNAMICJOBS = "TuneDynamicJobs"

const GRB_INT_PAR_UPDATEMODE = "UpdateMode"

const GRB_INT_PAR_OBJNUMBER = "ObjNumber"

const GRB_INT_PAR_MULTIOBJMETHOD = "MultiObjMethod"

const GRB_INT_PAR_MULTIOBJPRE = "MultiObjPre"

const GRB_INT_PAR_SCENARIONUMBER = "ScenarioNumber"

const GRB_INT_PAR_POOLSOLUTIONS = "PoolSolutions"

const GRB_DBL_PAR_POOLGAP = "PoolGap"

const GRB_DBL_PAR_POOLGAPABS = "PoolGapAbs"

const GRB_INT_PAR_POOLSEARCHMODE = "PoolSearchMode"

const GRB_INT_PAR_IGNORENAMES = "IgnoreNames"

const GRB_INT_PAR_STARTNUMBER = "StartNumber"

const GRB_INT_PAR_PARTITIONPLACE = "PartitionPlace"

const GRB_INT_PAR_FUNCPIECES = "FuncPieces"

const GRB_DBL_PAR_FUNCPIECELENGTH = "FuncPieceLength"

const GRB_DBL_PAR_FUNCPIECEERROR = "FuncPieceError"

const GRB_DBL_PAR_FUNCPIECERATIO = "FuncPieceRatio"

const GRB_DBL_PAR_FUNCMAXVAL = "FuncMaxVal"

const GRB_INT_PAR_FUNCNONLINEAR = "FuncNonlinear"

const GRB_STR_PAR_DUMMY = "Dummy"

const GRB_STR_PAR_JOBID = "JobID"

const GRB_CUTS_AUTO = -1

const GRB_CUTS_OFF = 0

const GRB_CUTS_CONSERVATIVE = 1

const GRB_CUTS_AGGRESSIVE = 2

const GRB_CUTS_VERYAGGRESSIVE = 3

const GRB_PRESOLVE_AUTO = -1

const GRB_PRESOLVE_OFF = 0

const GRB_PRESOLVE_CONSERVATIVE = 1

const GRB_PRESOLVE_AGGRESSIVE = 2

const GRB_METHOD_NONE = -1

const GRB_METHOD_AUTO = -1

const GRB_METHOD_PRIMAL = 0

const GRB_METHOD_DUAL = 1

const GRB_METHOD_BARRIER = 2

const GRB_METHOD_CONCURRENT = 3

const GRB_METHOD_DETERMINISTIC_CONCURRENT = 4

const GRB_METHOD_DETERMINISTIC_CONCURRENT_SIMPLEX = 5

const GRB_CONCURRENTMETHOD_AUTO = -1

const GRB_CONCURRENTMETHOD_BARRIER_PRIMAL_DUAL = 0

const GRB_CONCURRENTMETHOD_BARRIER_DUAL = 1

const GRB_CONCURRENTMETHOD_BARRIER_PRIMAL = 2

const GRB_CONCURRENTMETHOD_PRIMAL_DUAL = 3

const GRB_BARHOMOGENEOUS_AUTO = -1

const GRB_BARHOMOGENEOUS_OFF = 0

const GRB_BARHOMOGENEOUS_ON = 1

const GRB_BARORDER_AUTOMATIC = -1

const GRB_BARORDER_AMD = 0

const GRB_BARORDER_NESTEDDISSECTION = 1

const GRB_MIPFOCUS_BALANCED = 0

const GRB_MIPFOCUS_FEASIBILITY = 1

const GRB_MIPFOCUS_OPTIMALITY = 2

const GRB_MIPFOCUS_BESTBOUND = 3

const GRB_SIMPLEXPRICING_AUTO = -1

const GRB_SIMPLEXPRICING_PARTIAL = 0

const GRB_SIMPLEXPRICING_STEEPEST_EDGE = 1

const GRB_SIMPLEXPRICING_DEVEX = 2

const GRB_SIMPLEXPRICING_STEEPEST_QUICK = 3

const GRB_VARBRANCH_AUTO = -1

const GRB_VARBRANCH_PSEUDO_REDUCED = 0

const GRB_VARBRANCH_PSEUDO_SHADOW = 1

const GRB_VARBRANCH_MAX_INFEAS = 2

const GRB_VARBRANCH_STRONG = 3

const GRB_PARTITION_EARLY = 16

const GRB_PARTITION_ROOTSTART = 8

const GRB_PARTITION_ROOTEND = 4

const GRB_PARTITION_NODES = 2

const GRB_PARTITION_CLEANUP = 1

const GRB_PHASE_MIP_NOREL = 0

const GRB_PHASE_MIP_SEARCH = 1

const GRB_PHASE_MIP_IMPROVE = 2

function GRBloadenv(envP, logfilename)
    return GRBloadenvinternal(
        envP,
        logfilename,
        GRB_VERSION_MAJOR,
        GRB_VERSION_MINOR,
        GRB_VERSION_TECHNICAL,
    )
end

function GRBemptyenv(envP)
    return GRBemptyenvinternal(
        envP,
        GRB_VERSION_MAJOR,
        GRB_VERSION_MINOR,
        GRB_VERSION_TECHNICAL,
    )
end

const GRB_BATCH_STATUS_UNKNOWN = 0

const GRB_BATCH_CREATED = 1

const GRB_BATCH_SUBMITTED = 2

const GRB_BATCH_ABORTED = 3

const GRB_BATCH_FAILED = 4

const GRB_BATCH_COMPLETED = 5
