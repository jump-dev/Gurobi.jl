# Julia wrapper for header: gurobi_c.h
# Automatically generated using Clang.jl


function GRBgetattrinfo(model, attrname, datatypeP, attrtypeP, settableP)
    ccall((:GRBgetattrinfo, libgurobi), Cint, (Ptr{GRBmodel}, Ptr{Cchar}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}), model, attrname, datatypeP, attrtypeP, settableP)
end

function GRBisattravailable(model, attrname)
    ccall((:GRBisattravailable, libgurobi), Cint, (Ptr{GRBmodel}, Ptr{Cchar}), model, attrname)
end

function GRBgetintattr(model, attrname, valueP)
    ccall((:GRBgetintattr, libgurobi), Cint, (Ptr{GRBmodel}, Ptr{Cchar}, Ptr{Cint}), model, attrname, valueP)
end

function GRBsetintattr(model, attrname, newvalue)
    ccall((:GRBsetintattr, libgurobi), Cint, (Ptr{GRBmodel}, Ptr{Cchar}, Cint), model, attrname, newvalue)
end

function GRBgetintattrelement(model, attrname, element, valueP)
    ccall((:GRBgetintattrelement, libgurobi), Cint, (Ptr{GRBmodel}, Ptr{Cchar}, Cint, Ptr{Cint}), model, attrname, element, valueP)
end

function GRBsetintattrelement(model, attrname, element, newvalue)
    ccall((:GRBsetintattrelement, libgurobi), Cint, (Ptr{GRBmodel}, Ptr{Cchar}, Cint, Cint), model, attrname, element, newvalue)
end

function GRBgetintattrarray(model, attrname, first, len, values)
    ccall((:GRBgetintattrarray, libgurobi), Cint, (Ptr{GRBmodel}, Ptr{Cchar}, Cint, Cint, Ptr{Cint}), model, attrname, first, len, values)
end

function GRBsetintattrarray(model, attrname, first, len, newvalues)
    ccall((:GRBsetintattrarray, libgurobi), Cint, (Ptr{GRBmodel}, Ptr{Cchar}, Cint, Cint, Ptr{Cint}), model, attrname, first, len, newvalues)
end

function GRBgetintattrlist(model, attrname, len, ind, values)
    ccall((:GRBgetintattrlist, libgurobi), Cint, (Ptr{GRBmodel}, Ptr{Cchar}, Cint, Ptr{Cint}, Ptr{Cint}), model, attrname, len, ind, values)
end

function GRBsetintattrlist(model, attrname, len, ind, newvalues)
    ccall((:GRBsetintattrlist, libgurobi), Cint, (Ptr{GRBmodel}, Ptr{Cchar}, Cint, Ptr{Cint}, Ptr{Cint}), model, attrname, len, ind, newvalues)
end

function GRBgetcharattrelement(model, attrname, element, valueP)
    ccall((:GRBgetcharattrelement, libgurobi), Cint, (Ptr{GRBmodel}, Ptr{Cchar}, Cint, Ptr{Cchar}), model, attrname, element, valueP)
end

function GRBsetcharattrelement(model, attrname, element, newvalue)
    ccall((:GRBsetcharattrelement, libgurobi), Cint, (Ptr{GRBmodel}, Ptr{Cchar}, Cint, UInt8), model, attrname, element, newvalue)
end

function GRBgetcharattrarray(model, attrname, first, len, values)
    ccall((:GRBgetcharattrarray, libgurobi), Cint, (Ptr{GRBmodel}, Ptr{Cchar}, Cint, Cint, Ptr{Cchar}), model, attrname, first, len, values)
end

function GRBsetcharattrarray(model, attrname, first, len, newvalues)
    ccall((:GRBsetcharattrarray, libgurobi), Cint, (Ptr{GRBmodel}, Ptr{Cchar}, Cint, Cint, Ptr{Cchar}), model, attrname, first, len, newvalues)
end

function GRBgetcharattrlist(model, attrname, len, ind, values)
    ccall((:GRBgetcharattrlist, libgurobi), Cint, (Ptr{GRBmodel}, Ptr{Cchar}, Cint, Ptr{Cint}, Ptr{Cchar}), model, attrname, len, ind, values)
end

function GRBsetcharattrlist(model, attrname, len, ind, newvalues)
    ccall((:GRBsetcharattrlist, libgurobi), Cint, (Ptr{GRBmodel}, Ptr{Cchar}, Cint, Ptr{Cint}, Ptr{Cchar}), model, attrname, len, ind, newvalues)
end

function GRBgetdblattr(model, attrname, valueP)
    ccall((:GRBgetdblattr, libgurobi), Cint, (Ptr{GRBmodel}, Ptr{Cchar}, Ptr{Cdouble}), model, attrname, valueP)
end

function GRBsetdblattr(model, attrname, newvalue)
    ccall((:GRBsetdblattr, libgurobi), Cint, (Ptr{GRBmodel}, Ptr{Cchar}, Cdouble), model, attrname, newvalue)
end

function GRBgetdblattrelement(model, attrname, element, valueP)
    ccall((:GRBgetdblattrelement, libgurobi), Cint, (Ptr{GRBmodel}, Ptr{Cchar}, Cint, Ptr{Cdouble}), model, attrname, element, valueP)
end

function GRBsetdblattrelement(model, attrname, element, newvalue)
    ccall((:GRBsetdblattrelement, libgurobi), Cint, (Ptr{GRBmodel}, Ptr{Cchar}, Cint, Cdouble), model, attrname, element, newvalue)
end

function GRBgetdblattrarray(model, attrname, first, len, values)
    ccall((:GRBgetdblattrarray, libgurobi), Cint, (Ptr{GRBmodel}, Ptr{Cchar}, Cint, Cint, Ptr{Cdouble}), model, attrname, first, len, values)
end

function GRBsetdblattrarray(model, attrname, first, len, newvalues)
    ccall((:GRBsetdblattrarray, libgurobi), Cint, (Ptr{GRBmodel}, Ptr{Cchar}, Cint, Cint, Ptr{Cdouble}), model, attrname, first, len, newvalues)
end

function GRBgetdblattrlist(model, attrname, len, ind, values)
    ccall((:GRBgetdblattrlist, libgurobi), Cint, (Ptr{GRBmodel}, Ptr{Cchar}, Cint, Ptr{Cint}, Ptr{Cdouble}), model, attrname, len, ind, values)
end

function GRBsetdblattrlist(model, attrname, len, ind, newvalues)
    ccall((:GRBsetdblattrlist, libgurobi), Cint, (Ptr{GRBmodel}, Ptr{Cchar}, Cint, Ptr{Cint}, Ptr{Cdouble}), model, attrname, len, ind, newvalues)
end

function GRBgetstrattr(model, attrname, valueP)
    ccall((:GRBgetstrattr, libgurobi), Cint, (Ptr{GRBmodel}, Ptr{Cchar}, Ptr{Ptr{Cchar}}), model, attrname, valueP)
end

function GRBsetstrattr(model, attrname, newvalue)
    ccall((:GRBsetstrattr, libgurobi), Cint, (Ptr{GRBmodel}, Ptr{Cchar}, Ptr{Cchar}), model, attrname, newvalue)
end

function GRBgetstrattrelement(model, attrname, element, valueP)
    ccall((:GRBgetstrattrelement, libgurobi), Cint, (Ptr{GRBmodel}, Ptr{Cchar}, Cint, Ptr{Ptr{Cchar}}), model, attrname, element, valueP)
end

function GRBsetstrattrelement(model, attrname, element, newvalue)
    ccall((:GRBsetstrattrelement, libgurobi), Cint, (Ptr{GRBmodel}, Ptr{Cchar}, Cint, Ptr{Cchar}), model, attrname, element, newvalue)
end

function GRBgetstrattrarray(model, attrname, first, len, values)
    ccall((:GRBgetstrattrarray, libgurobi), Cint, (Ptr{GRBmodel}, Ptr{Cchar}, Cint, Cint, Ptr{Ptr{Cchar}}), model, attrname, first, len, values)
end

function GRBsetstrattrarray(model, attrname, first, len, newvalues)
    ccall((:GRBsetstrattrarray, libgurobi), Cint, (Ptr{GRBmodel}, Ptr{Cchar}, Cint, Cint, Ptr{Ptr{Cchar}}), model, attrname, first, len, newvalues)
end

function GRBgetstrattrlist(model, attrname, len, ind, values)
    ccall((:GRBgetstrattrlist, libgurobi), Cint, (Ptr{GRBmodel}, Ptr{Cchar}, Cint, Ptr{Cint}, Ptr{Ptr{Cchar}}), model, attrname, len, ind, values)
end

function GRBsetstrattrlist(model, attrname, len, ind, newvalues)
    ccall((:GRBsetstrattrlist, libgurobi), Cint, (Ptr{GRBmodel}, Ptr{Cchar}, Cint, Ptr{Cint}, Ptr{Ptr{Cchar}}), model, attrname, len, ind, newvalues)
end

function GRBsetcallbackfunc(model, cb, usrdata)
    ccall((:GRBsetcallbackfunc, libgurobi), Cint, (Ptr{GRBmodel}, Ptr{Cvoid}, Ptr{Cvoid}), model, cb, usrdata)
end

function GRBgetcallbackfunc(model, cbP)
    ccall((:GRBgetcallbackfunc, libgurobi), Cint, (Ptr{GRBmodel}, Ptr{Ptr{Cvoid}}), model, cbP)
end

function GRBsetlogcallbackfunc(model, cb)
    ccall((:GRBsetlogcallbackfunc, libgurobi), Cint, (Ptr{GRBmodel}, Ptr{Cvoid}), model, cb)
end

function GRBsetlogcallbackfuncenv(env, cb)
    ccall((:GRBsetlogcallbackfuncenv, libgurobi), Cint, (Ptr{GRBenv}, Ptr{Cvoid}), env, cb)
end

function GRBcbget(cbdata, where, what, resultP)
    ccall((:GRBcbget, libgurobi), Cint, (Ptr{Cvoid}, Cint, Cint, Ptr{Cvoid}), cbdata, where, what, resultP)
end

function GRBcbsetparam(cbdata, paramname, newvalue)
    ccall((:GRBcbsetparam, libgurobi), Cint, (Ptr{Cvoid}, Ptr{Cchar}, Ptr{Cchar}), cbdata, paramname, newvalue)
end

function GRBcbsolution(cbdata, solution, objvalP)
    ccall((:GRBcbsolution, libgurobi), Cint, (Ptr{Cvoid}, Ptr{Cdouble}, Ptr{Cdouble}), cbdata, solution, objvalP)
end

function GRBcbcut(cbdata, cutlen, cutind, cutval, cutsense, cutrhs)
    ccall((:GRBcbcut, libgurobi), Cint, (Ptr{Cvoid}, Cint, Ptr{Cint}, Ptr{Cdouble}, UInt8, Cdouble), cbdata, cutlen, cutind, cutval, cutsense, cutrhs)
end

function GRBcblazy(cbdata, lazylen, lazyind, lazyval, lazysense, lazyrhs)
    ccall((:GRBcblazy, libgurobi), Cint, (Ptr{Cvoid}, Cint, Ptr{Cint}, Ptr{Cdouble}, UInt8, Cdouble), cbdata, lazylen, lazyind, lazyval, lazysense, lazyrhs)
end

function GRBgetcoeff(model, constr, var, valP)
    ccall((:GRBgetcoeff, libgurobi), Cint, (Ptr{GRBmodel}, Cint, Cint, Ptr{Cdouble}), model, constr, var, valP)
end

function GRBgetconstrs(model, numnzP, cbeg, cind, cval, start, len)
    ccall((:GRBgetconstrs, libgurobi), Cint, (Ptr{GRBmodel}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}, Cint, Cint), model, numnzP, cbeg, cind, cval, start, len)
end

function GRBXgetconstrs(model, numnzP, cbeg, cind, cval, start, len)
    ccall((:GRBXgetconstrs, libgurobi), Cint, (Ptr{GRBmodel}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}, Cint, Cint), model, numnzP, cbeg, cind, cval, start, len)
end

function GRBgetvars(model, numnzP, vbeg, vind, vval, start, len)
    ccall((:GRBgetvars, libgurobi), Cint, (Ptr{GRBmodel}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}, Cint, Cint), model, numnzP, vbeg, vind, vval, start, len)
end

function GRBXgetvars(model, numnzP, vbeg, vind, vval, start, len)
    ccall((:GRBXgetvars, libgurobi), Cint, (Ptr{GRBmodel}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}, Cint, Cint), model, numnzP, vbeg, vind, vval, start, len)
end

function GRBgetsos(model, nummembersP, sostype, beg, ind, weight, start, len)
    ccall((:GRBgetsos, libgurobi), Cint, (Ptr{GRBmodel}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}, Cint, Cint), model, nummembersP, sostype, beg, ind, weight, start, len)
end

function GRBgetgenconstrMax(model, genconstr, resvarP, nvarsP, vars, constantP)
    ccall((:GRBgetgenconstrMax, libgurobi), Cint, (Ptr{GRBmodel}, Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}), model, genconstr, resvarP, nvarsP, vars, constantP)
end

function GRBgetgenconstrMin(model, genconstr, resvarP, nvarsP, vars, constantP)
    ccall((:GRBgetgenconstrMin, libgurobi), Cint, (Ptr{GRBmodel}, Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}), model, genconstr, resvarP, nvarsP, vars, constantP)
end

function GRBgetgenconstrAbs(model, genconstr, resvarP, argvarP)
    ccall((:GRBgetgenconstrAbs, libgurobi), Cint, (Ptr{GRBmodel}, Cint, Ptr{Cint}, Ptr{Cint}), model, genconstr, resvarP, argvarP)
end

function GRBgetgenconstrAnd(model, genconstr, resvarP, nvarsP, vars)
    ccall((:GRBgetgenconstrAnd, libgurobi), Cint, (Ptr{GRBmodel}, Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}), model, genconstr, resvarP, nvarsP, vars)
end

function GRBgetgenconstrOr(model, genconstr, resvarP, nvarsP, vars)
    ccall((:GRBgetgenconstrOr, libgurobi), Cint, (Ptr{GRBmodel}, Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}), model, genconstr, resvarP, nvarsP, vars)
end

function GRBgetgenconstrIndicator(model, genconstr, binvarP, binvalP, nvarsP, vars, vals, senseP, rhsP)
    ccall((:GRBgetgenconstrIndicator, libgurobi), Cint, (Ptr{GRBmodel}, Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cchar}, Ptr{Cdouble}), model, genconstr, binvarP, binvalP, nvarsP, vars, vals, senseP, rhsP)
end

function GRBgetgenconstrPWL(model, genconstr, xvarP, yvarP, nptsP, xpts, ypts)
    ccall((:GRBgetgenconstrPWL, libgurobi), Cint, (Ptr{GRBmodel}, Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cdouble}), model, genconstr, xvarP, yvarP, nptsP, xpts, ypts)
end

function GRBgetgenconstrPoly(model, genconstr, xvarP, yvarP, plenP, p)
    ccall((:GRBgetgenconstrPoly, libgurobi), Cint, (Ptr{GRBmodel}, Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}), model, genconstr, xvarP, yvarP, plenP, p)
end

function GRBgetgenconstrExpA(model, genconstr, xvarP, yvarP, aP)
    ccall((:GRBgetgenconstrExpA, libgurobi), Cint, (Ptr{GRBmodel}, Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}), model, genconstr, xvarP, yvarP, aP)
end

function GRBgetgenconstrLogA(model, genconstr, xvarP, yvarP, aP)
    ccall((:GRBgetgenconstrLogA, libgurobi), Cint, (Ptr{GRBmodel}, Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}), model, genconstr, xvarP, yvarP, aP)
end

function GRBgetgenconstrPow(model, genconstr, xvarP, yvarP, aP)
    ccall((:GRBgetgenconstrPow, libgurobi), Cint, (Ptr{GRBmodel}, Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}), model, genconstr, xvarP, yvarP, aP)
end

function GRBgetgenconstrExp(model, genconstr, xvarP, yvarP)
    ccall((:GRBgetgenconstrExp, libgurobi), Cint, (Ptr{GRBmodel}, Cint, Ptr{Cint}, Ptr{Cint}), model, genconstr, xvarP, yvarP)
end

function GRBgetgenconstrLog(model, genconstr, xvarP, yvarP)
    ccall((:GRBgetgenconstrLog, libgurobi), Cint, (Ptr{GRBmodel}, Cint, Ptr{Cint}, Ptr{Cint}), model, genconstr, xvarP, yvarP)
end

function GRBgetgenconstrSin(model, genconstr, xvarP, yvarP)
    ccall((:GRBgetgenconstrSin, libgurobi), Cint, (Ptr{GRBmodel}, Cint, Ptr{Cint}, Ptr{Cint}), model, genconstr, xvarP, yvarP)
end

function GRBgetgenconstrCos(model, genconstr, xvarP, yvarP)
    ccall((:GRBgetgenconstrCos, libgurobi), Cint, (Ptr{GRBmodel}, Cint, Ptr{Cint}, Ptr{Cint}), model, genconstr, xvarP, yvarP)
end

function GRBgetgenconstrTan(model, genconstr, xvarP, yvarP)
    ccall((:GRBgetgenconstrTan, libgurobi), Cint, (Ptr{GRBmodel}, Cint, Ptr{Cint}, Ptr{Cint}), model, genconstr, xvarP, yvarP)
end

function GRBgetq(model, numqnzP, qrow, qcol, qval)
    ccall((:GRBgetq, libgurobi), Cint, (Ptr{GRBmodel}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}), model, numqnzP, qrow, qcol, qval)
end

function GRBgetqconstr(model, qconstr, numlnzP, lind, lval, numqnzP, qrow, qcol, qval)
    ccall((:GRBgetqconstr, libgurobi), Cint, (Ptr{GRBmodel}, Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}), model, qconstr, numlnzP, lind, lval, numqnzP, qrow, qcol, qval)
end

function GRBgetvarbyname(model, name, indexP)
    ccall((:GRBgetvarbyname, libgurobi), Cint, (Ptr{GRBmodel}, Ptr{Cchar}, Ptr{Cint}), model, name, indexP)
end

function GRBgetconstrbyname(model, name, indexP)
    ccall((:GRBgetconstrbyname, libgurobi), Cint, (Ptr{GRBmodel}, Ptr{Cchar}, Ptr{Cint}), model, name, indexP)
end

function GRBgetqconstrbyname(model, name, indexP)
    ccall((:GRBgetqconstrbyname, libgurobi), Cint, (Ptr{GRBmodel}, Ptr{Cchar}, Ptr{Cint}), model, name, indexP)
end

function GRBgetpwlobj(model, var, pointsP, x, y)
    ccall((:GRBgetpwlobj, libgurobi), Cint, (Ptr{GRBmodel}, Cint, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cdouble}), model, var, pointsP, x, y)
end

function GRBoptimize(model)
    ccall((:GRBoptimize, libgurobi), Cint, (Ptr{GRBmodel},), model)
end

function GRBoptimizeasync(model)
    ccall((:GRBoptimizeasync, libgurobi), Cint, (Ptr{GRBmodel},), model)
end

function GRBoptimizebatch(model, batchid)
    ccall((:GRBoptimizebatch, libgurobi), Cint, (Ptr{GRBmodel}, Ptr{Cchar}), model, batchid)
end

function GRBcopymodel(model)
    ccall((:GRBcopymodel, libgurobi), Ptr{GRBmodel}, (Ptr{GRBmodel},), model)
end

function GRBfixmodel(model, fixedP)
    ccall((:GRBfixmodel, libgurobi), Cint, (Ptr{GRBmodel}, Ptr{Ptr{GRBmodel}}), model, fixedP)
end

function GRBfeasrelax(model, relaxobjtype, minrelax, lbpen, ubpen, rhspen, feasobjP)
    ccall((:GRBfeasrelax, libgurobi), Cint, (Ptr{GRBmodel}, Cint, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}), model, relaxobjtype, minrelax, lbpen, ubpen, rhspen, feasobjP)
end

function GRBsinglescenariomodel(model, singlescenarioP)
    ccall((:GRBsinglescenariomodel, libgurobi), Cint, (Ptr{GRBmodel}, Ptr{Ptr{GRBmodel}}), model, singlescenarioP)
end

function GRBgetcbwhatinfo(cbdata, what, typeP, sizeP)
    ccall((:GRBgetcbwhatinfo, libgurobi), Cint, (Ptr{Cvoid}, Cint, Ptr{Cint}, Ptr{Cint}), cbdata, what, typeP, sizeP)
end

function GRBrelaxmodel(model, relaxedP)
    ccall((:GRBrelaxmodel, libgurobi), Cint, (Ptr{GRBmodel}, Ptr{Ptr{GRBmodel}}), model, relaxedP)
end

function GRBconverttofixed(model)
    ccall((:GRBconverttofixed, libgurobi), Cint, (Ptr{GRBmodel},), model)
end

function GRBpresolvemodel(model, presolvedP)
    ccall((:GRBpresolvemodel, libgurobi), Cint, (Ptr{GRBmodel}, Ptr{Ptr{GRBmodel}}), model, presolvedP)
end

function GRBiismodel(model, iisP)
    ccall((:GRBiismodel, libgurobi), Cint, (Ptr{GRBmodel}, Ptr{Ptr{GRBmodel}}), model, iisP)
end

function GRBfeasibility(model, feasP)
    ccall((:GRBfeasibility, libgurobi), Cint, (Ptr{GRBmodel}, Ptr{Ptr{GRBmodel}}), model, feasP)
end

function GRBlinearizemodel(model, linearizedP)
    ccall((:GRBlinearizemodel, libgurobi), Cint, (Ptr{GRBmodel}, Ptr{Ptr{GRBmodel}}), model, linearizedP)
end

function GRBresultmodel(model, type, resultP)
    ccall((:GRBresultmodel, libgurobi), Cint, (Ptr{GRBmodel}, Ptr{Cchar}, Ptr{Ptr{GRBmodel}}), model, type, resultP)
end

function GRBfixedmodel(model)
    ccall((:GRBfixedmodel, libgurobi), Ptr{GRBmodel}, (Ptr{GRBmodel},), model)
end

function GRBloadenvsyscb(envP, logfilename, malloccb, calloccb, realloccb, freecb, threadcreatecb, threadjoincb, syscbusrdata)
    ccall((:GRBloadenvsyscb, libgurobi), Cint, (Ptr{Ptr{GRBenv}}, Ptr{Cchar}, Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}), envP, logfilename, malloccb, calloccb, realloccb, freecb, threadcreatecb, threadjoincb, syscbusrdata)
end

function GRBemptyenvadv(envP, apitype, major, minor, tech, malloccb, calloccb, realloccb, freecb, threadcreatecb, threadjoincb, syscbusrdata)
    ccall((:GRBemptyenvadv, libgurobi), Cint, (Ptr{Ptr{GRBenv}}, Cint, Cint, Cint, Cint, Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}), envP, apitype, major, minor, tech, malloccb, calloccb, realloccb, freecb, threadcreatecb, threadjoincb, syscbusrdata)
end

function GRBreadmodel(env, filename, modelP)
    ccall((:GRBreadmodel, libgurobi), Cint, (Ptr{GRBenv}, Ptr{Cchar}, Ptr{Ptr{GRBmodel}}), env, filename, modelP)
end

function GRBread(model, filename)
    ccall((:GRBread, libgurobi), Cint, (Ptr{GRBmodel}, Ptr{Cchar}), model, filename)
end

function GRBwrite(model, filename)
    ccall((:GRBwrite, libgurobi), Cint, (Ptr{GRBmodel}, Ptr{Cchar}), model, filename)
end

function GRBismodelfile(filename)
    ccall((:GRBismodelfile, libgurobi), Cint, (Ptr{Cchar},), filename)
end

function GRBisattrfile(filename)
    ccall((:GRBisattrfile, libgurobi), Cint, (Ptr{Cchar},), filename)
end

function GRBfiletype(filename)
    ccall((:GRBfiletype, libgurobi), Cint, (Ptr{Cchar},), filename)
end

function GRBisrecordfile(filename)
    ccall((:GRBisrecordfile, libgurobi), Cint, (Ptr{Cchar},), filename)
end

function GRBgetjsonsolution(model, buffP)
    ccall((:GRBgetjsonsolution, libgurobi), Cint, (Ptr{GRBmodel}, Ptr{Ptr{Cchar}}), model, buffP)
end

function GRBloadjson(env, fname, buffP)
    ccall((:GRBloadjson, libgurobi), Cint, (Ptr{GRBenv}, Ptr{Cchar}, Ptr{Ptr{Cchar}}), env, fname, buffP)
end

function GRBnewmodel(env, modelP, Pname, numvars, obj, lb, ub, vtype, varnames)
    ccall((:GRBnewmodel, libgurobi), Cint, (Ptr{GRBenv}, Ptr{Ptr{GRBmodel}}, Ptr{Cchar}, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cchar}, Ptr{Ptr{Cchar}}), env, modelP, Pname, numvars, obj, lb, ub, vtype, varnames)
end

function GRBloadmodel(env, modelP, Pname, numvars, numconstrs, objsense, objcon, obj, sense, rhs, vbeg, vlen, vind, vval, lb, ub, vtype, varnames, constrnames)
    ccall((:GRBloadmodel, libgurobi), Cint, (Ptr{GRBenv}, Ptr{Ptr{GRBmodel}}, Ptr{Cchar}, Cint, Cint, Cint, Cdouble, Ptr{Cdouble}, Ptr{Cchar}, Ptr{Cdouble}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cchar}, Ptr{Ptr{Cchar}}, Ptr{Ptr{Cchar}}), env, modelP, Pname, numvars, numconstrs, objsense, objcon, obj, sense, rhs, vbeg, vlen, vind, vval, lb, ub, vtype, varnames, constrnames)
end

function GRBXloadmodel(env, modelP, Pname, numvars, numconstrs, objsense, objcon, obj, sense, rhs, vbeg, vlen, vind, vval, lb, ub, vtype, varnames, constrnames)
    ccall((:GRBXloadmodel, libgurobi), Cint, (Ptr{GRBenv}, Ptr{Ptr{GRBmodel}}, Ptr{Cchar}, Cint, Cint, Cint, Cdouble, Ptr{Cdouble}, Ptr{Cchar}, Ptr{Cdouble}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cchar}, Ptr{Ptr{Cchar}}, Ptr{Ptr{Cchar}}), env, modelP, Pname, numvars, numconstrs, objsense, objcon, obj, sense, rhs, vbeg, vlen, vind, vval, lb, ub, vtype, varnames, constrnames)
end

function GRBaddvar(model, numnz, vind, vval, obj, lb, ub, vtype, varname)
    ccall((:GRBaddvar, libgurobi), Cint, (Ptr{GRBmodel}, Cint, Ptr{Cint}, Ptr{Cdouble}, Cdouble, Cdouble, Cdouble, UInt8, Ptr{Cchar}), model, numnz, vind, vval, obj, lb, ub, vtype, varname)
end

function GRBaddvars(model, numvars, numnz, vbeg, vind, vval, obj, lb, ub, vtype, varnames)
    ccall((:GRBaddvars, libgurobi), Cint, (Ptr{GRBmodel}, Cint, Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cchar}, Ptr{Ptr{Cchar}}), model, numvars, numnz, vbeg, vind, vval, obj, lb, ub, vtype, varnames)
end

function GRBXaddvars(model, numvars, numnz, vbeg, vind, vval, obj, lb, ub, vtype, varnames)
    ccall((:GRBXaddvars, libgurobi), Cint, (Ptr{GRBmodel}, Cint, Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cchar}, Ptr{Ptr{Cchar}}), model, numvars, numnz, vbeg, vind, vval, obj, lb, ub, vtype, varnames)
end

function GRBaddconstr(model, numnz, cind, cval, sense, rhs, constrname)
    ccall((:GRBaddconstr, libgurobi), Cint, (Ptr{GRBmodel}, Cint, Ptr{Cint}, Ptr{Cdouble}, UInt8, Cdouble, Ptr{Cchar}), model, numnz, cind, cval, sense, rhs, constrname)
end

function GRBaddconstrs(model, numconstrs, numnz, cbeg, cind, cval, sense, rhs, constrnames)
    ccall((:GRBaddconstrs, libgurobi), Cint, (Ptr{GRBmodel}, Cint, Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cchar}, Ptr{Cdouble}, Ptr{Ptr{Cchar}}), model, numconstrs, numnz, cbeg, cind, cval, sense, rhs, constrnames)
end

function GRBXaddconstrs(model, numconstrs, numnz, cbeg, cind, cval, sense, rhs, constrnames)
    ccall((:GRBXaddconstrs, libgurobi), Cint, (Ptr{GRBmodel}, Cint, Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cchar}, Ptr{Cdouble}, Ptr{Ptr{Cchar}}), model, numconstrs, numnz, cbeg, cind, cval, sense, rhs, constrnames)
end

function GRBaddrangeconstr(model, numnz, cind, cval, lower, upper, constrname)
    ccall((:GRBaddrangeconstr, libgurobi), Cint, (Ptr{GRBmodel}, Cint, Ptr{Cint}, Ptr{Cdouble}, Cdouble, Cdouble, Ptr{Cchar}), model, numnz, cind, cval, lower, upper, constrname)
end

function GRBaddrangeconstrs(model, numconstrs, numnz, cbeg, cind, cval, lower, upper, constrnames)
    ccall((:GRBaddrangeconstrs, libgurobi), Cint, (Ptr{GRBmodel}, Cint, Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Ptr{Cchar}}), model, numconstrs, numnz, cbeg, cind, cval, lower, upper, constrnames)
end

function GRBXaddrangeconstrs(model, numconstrs, numnz, cbeg, cind, cval, lower, upper, constrnames)
    ccall((:GRBXaddrangeconstrs, libgurobi), Cint, (Ptr{GRBmodel}, Cint, Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Ptr{Cchar}}), model, numconstrs, numnz, cbeg, cind, cval, lower, upper, constrnames)
end

function GRBaddsos(model, numsos, nummembers, types, beg, ind, weight)
    ccall((:GRBaddsos, libgurobi), Cint, (Ptr{GRBmodel}, Cint, Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}), model, numsos, nummembers, types, beg, ind, weight)
end

function GRBaddgenconstrMax(model, name, resvar, nvars, vars, constant)
    ccall((:GRBaddgenconstrMax, libgurobi), Cint, (Ptr{GRBmodel}, Ptr{Cchar}, Cint, Cint, Ptr{Cint}, Cdouble), model, name, resvar, nvars, vars, constant)
end

function GRBaddgenconstrMin(model, name, resvar, nvars, vars, constant)
    ccall((:GRBaddgenconstrMin, libgurobi), Cint, (Ptr{GRBmodel}, Ptr{Cchar}, Cint, Cint, Ptr{Cint}, Cdouble), model, name, resvar, nvars, vars, constant)
end

function GRBaddgenconstrAbs(model, name, resvar, argvar)
    ccall((:GRBaddgenconstrAbs, libgurobi), Cint, (Ptr{GRBmodel}, Ptr{Cchar}, Cint, Cint), model, name, resvar, argvar)
end

function GRBaddgenconstrAnd(model, name, resvar, nvars, vars)
    ccall((:GRBaddgenconstrAnd, libgurobi), Cint, (Ptr{GRBmodel}, Ptr{Cchar}, Cint, Cint, Ptr{Cint}), model, name, resvar, nvars, vars)
end

function GRBaddgenconstrOr(model, name, resvar, nvars, vars)
    ccall((:GRBaddgenconstrOr, libgurobi), Cint, (Ptr{GRBmodel}, Ptr{Cchar}, Cint, Cint, Ptr{Cint}), model, name, resvar, nvars, vars)
end

function GRBaddgenconstrIndicator(model, name, binvar, binval, nvars, vars, vals, sense, rhs)
    ccall((:GRBaddgenconstrIndicator, libgurobi), Cint, (Ptr{GRBmodel}, Ptr{Cchar}, Cint, Cint, Cint, Ptr{Cint}, Ptr{Cdouble}, UInt8, Cdouble), model, name, binvar, binval, nvars, vars, vals, sense, rhs)
end

function GRBaddgenconstrPWL(model, name, xvar, yvar, npts, xpts, ypts)
    ccall((:GRBaddgenconstrPWL, libgurobi), Cint, (Ptr{GRBmodel}, Ptr{Cchar}, Cint, Cint, Cint, Ptr{Cdouble}, Ptr{Cdouble}), model, name, xvar, yvar, npts, xpts, ypts)
end

function GRBaddgenconstrPoly(model, name, xvar, yvar, plen, p, options)
    ccall((:GRBaddgenconstrPoly, libgurobi), Cint, (Ptr{GRBmodel}, Ptr{Cchar}, Cint, Cint, Cint, Ptr{Cdouble}, Ptr{Cchar}), model, name, xvar, yvar, plen, p, options)
end

function GRBaddgenconstrExpA(model, name, xvar, yvar, a, options)
    ccall((:GRBaddgenconstrExpA, libgurobi), Cint, (Ptr{GRBmodel}, Ptr{Cchar}, Cint, Cint, Cdouble, Ptr{Cchar}), model, name, xvar, yvar, a, options)
end

function GRBaddgenconstrLogA(model, name, xvar, yvar, a, options)
    ccall((:GRBaddgenconstrLogA, libgurobi), Cint, (Ptr{GRBmodel}, Ptr{Cchar}, Cint, Cint, Cdouble, Ptr{Cchar}), model, name, xvar, yvar, a, options)
end

function GRBaddgenconstrPow(model, name, xvar, yvar, a, options)
    ccall((:GRBaddgenconstrPow, libgurobi), Cint, (Ptr{GRBmodel}, Ptr{Cchar}, Cint, Cint, Cdouble, Ptr{Cchar}), model, name, xvar, yvar, a, options)
end

function GRBaddgenconstrExp(model, name, xvar, yvar, options)
    ccall((:GRBaddgenconstrExp, libgurobi), Cint, (Ptr{GRBmodel}, Ptr{Cchar}, Cint, Cint, Ptr{Cchar}), model, name, xvar, yvar, options)
end

function GRBaddgenconstrLog(model, name, xvar, yvar, options)
    ccall((:GRBaddgenconstrLog, libgurobi), Cint, (Ptr{GRBmodel}, Ptr{Cchar}, Cint, Cint, Ptr{Cchar}), model, name, xvar, yvar, options)
end

function GRBaddgenconstrSin(model, name, xvar, yvar, options)
    ccall((:GRBaddgenconstrSin, libgurobi), Cint, (Ptr{GRBmodel}, Ptr{Cchar}, Cint, Cint, Ptr{Cchar}), model, name, xvar, yvar, options)
end

function GRBaddgenconstrCos(model, name, xvar, yvar, options)
    ccall((:GRBaddgenconstrCos, libgurobi), Cint, (Ptr{GRBmodel}, Ptr{Cchar}, Cint, Cint, Ptr{Cchar}), model, name, xvar, yvar, options)
end

function GRBaddgenconstrTan(model, name, xvar, yvar, options)
    ccall((:GRBaddgenconstrTan, libgurobi), Cint, (Ptr{GRBmodel}, Ptr{Cchar}, Cint, Cint, Ptr{Cchar}), model, name, xvar, yvar, options)
end

function GRBaddqconstr(model, numlnz, lind, lval, numqnz, qrow, qcol, qval, sense, rhs, QCname)
    ccall((:GRBaddqconstr, libgurobi), Cint, (Ptr{GRBmodel}, Cint, Ptr{Cint}, Ptr{Cdouble}, Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}, UInt8, Cdouble, Ptr{Cchar}), model, numlnz, lind, lval, numqnz, qrow, qcol, qval, sense, rhs, QCname)
end

function GRBaddcone(model, nummembers, members)
    ccall((:GRBaddcone, libgurobi), Cint, (Ptr{GRBmodel}, Cint, Ptr{Cint}), model, nummembers, members)
end

function GRBaddqpterms(model, numqnz, qrow, qcol, qval)
    ccall((:GRBaddqpterms, libgurobi), Cint, (Ptr{GRBmodel}, Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}), model, numqnz, qrow, qcol, qval)
end

function GRBdelvars(model, len, ind)
    ccall((:GRBdelvars, libgurobi), Cint, (Ptr{GRBmodel}, Cint, Ptr{Cint}), model, len, ind)
end

function GRBdelconstrs(model, len, ind)
    ccall((:GRBdelconstrs, libgurobi), Cint, (Ptr{GRBmodel}, Cint, Ptr{Cint}), model, len, ind)
end

function GRBdelsos(model, len, ind)
    ccall((:GRBdelsos, libgurobi), Cint, (Ptr{GRBmodel}, Cint, Ptr{Cint}), model, len, ind)
end

function GRBdelgenconstrs(model, len, ind)
    ccall((:GRBdelgenconstrs, libgurobi), Cint, (Ptr{GRBmodel}, Cint, Ptr{Cint}), model, len, ind)
end

function GRBdelqconstrs(model, len, ind)
    ccall((:GRBdelqconstrs, libgurobi), Cint, (Ptr{GRBmodel}, Cint, Ptr{Cint}), model, len, ind)
end

function GRBdelq(model)
    ccall((:GRBdelq, libgurobi), Cint, (Ptr{GRBmodel},), model)
end

function GRBchgcoeffs(model, cnt, cind, vind, val)
    ccall((:GRBchgcoeffs, libgurobi), Cint, (Ptr{GRBmodel}, Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}), model, cnt, cind, vind, val)
end

function GRBXchgcoeffs(model, cnt, cind, vind, val)
    ccall((:GRBXchgcoeffs, libgurobi), Cint, (Ptr{GRBmodel}, Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}), model, cnt, cind, vind, val)
end

function GRBsetpwlobj(model, var, points, x, y)
    ccall((:GRBsetpwlobj, libgurobi), Cint, (Ptr{GRBmodel}, Cint, Cint, Ptr{Cdouble}, Ptr{Cdouble}), model, var, points, x, y)
end

function GRBupdatemodel(model)
    ccall((:GRBupdatemodel, libgurobi), Cint, (Ptr{GRBmodel},), model)
end

function GRBreset(model, clearall)
    ccall((:GRBreset, libgurobi), Cint, (Ptr{GRBmodel}, Cint), model, clearall)
end

function GRBresetmodel(model)
    ccall((:GRBresetmodel, libgurobi), Cint, (Ptr{GRBmodel},), model)
end

function GRBfreemodel(model)
    ccall((:GRBfreemodel, libgurobi), Cint, (Ptr{GRBmodel},), model)
end

function GRBcomputeIIS(model)
    ccall((:GRBcomputeIIS, libgurobi), Cint, (Ptr{GRBmodel},), model)
end

function GRBFSolve(model, b, x)
    ccall((:GRBFSolve, libgurobi), Cint, (Ptr{GRBmodel}, Ptr{GRBsvec}, Ptr{GRBsvec}), model, b, x)
end

function GRBBinvColj(model, j, x)
    ccall((:GRBBinvColj, libgurobi), Cint, (Ptr{GRBmodel}, Cint, Ptr{GRBsvec}), model, j, x)
end

function GRBBinvj(model, j, x)
    ccall((:GRBBinvj, libgurobi), Cint, (Ptr{GRBmodel}, Cint, Ptr{GRBsvec}), model, j, x)
end

function GRBBSolve(model, b, x)
    ccall((:GRBBSolve, libgurobi), Cint, (Ptr{GRBmodel}, Ptr{GRBsvec}, Ptr{GRBsvec}), model, b, x)
end

function GRBBinvi(model, i, x)
    ccall((:GRBBinvi, libgurobi), Cint, (Ptr{GRBmodel}, Cint, Ptr{GRBsvec}), model, i, x)
end

function GRBBinvRowi(model, i, x)
    ccall((:GRBBinvRowi, libgurobi), Cint, (Ptr{GRBmodel}, Cint, Ptr{GRBsvec}), model, i, x)
end

function GRBgetBasisHead(model, bhead)
    ccall((:GRBgetBasisHead, libgurobi), Cint, (Ptr{GRBmodel}, Ptr{Cint}), model, bhead)
end

function GRBcbstoponemultiobj(model, cbdata, objnum)
    ccall((:GRBcbstoponemultiobj, libgurobi), Cint, (Ptr{GRBmodel}, Ptr{Cvoid}, Cint), model, cbdata, objnum)
end

function GRBstrongbranch(model, num, cand, downobjbd, upobjbd, statusP)
    ccall((:GRBstrongbranch, libgurobi), Cint, (Ptr{GRBmodel}, Cint, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint}), model, num, cand, downobjbd, upobjbd, statusP)
end

function GRBcheckmodel(model)
    ccall((:GRBcheckmodel, libgurobi), Cint, (Ptr{GRBmodel},), model)
end

function GRBsetsignal(model)
    ccall((:GRBsetsignal, libgurobi), Cvoid, (Ptr{GRBmodel},), model)
end

function GRBterminate(model)
    ccall((:GRBterminate, libgurobi), Cvoid, (Ptr{GRBmodel},), model)
end

function GRBreplay(filename)
    ccall((:GRBreplay, libgurobi), Cint, (Ptr{Cchar},), filename)
end

function GRBsetobjective(model, sense, constant, lnz, lind, lval, qnz, qrow, qcol, qval)
    ccall((:GRBsetobjective, libgurobi), Cint, (Ptr{GRBmodel}, Cint, Cdouble, Cint, Ptr{Cint}, Ptr{Cdouble}, Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}), model, sense, constant, lnz, lind, lval, qnz, qrow, qcol, qval)
end

function GRBsetobjectiven(model, index, priority, weight, abstol, reltol, name, constant, lnz, lind, lval)
    ccall((:GRBsetobjectiven, libgurobi), Cint, (Ptr{GRBmodel}, Cint, Cint, Cdouble, Cdouble, Cdouble, Ptr{Cchar}, Cdouble, Cint, Ptr{Cint}, Ptr{Cdouble}), model, index, priority, weight, abstol, reltol, name, constant, lnz, lind, lval)
end

function GRBclean2(lenP, ind, val)
    ccall((:GRBclean2, libgurobi), Cvoid, (Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}), lenP, ind, val)
end

function GRBclean3(lenP, ind0, ind1, val)
    ccall((:GRBclean3, libgurobi), Cvoid, (Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}), lenP, ind0, ind1, val)
end

function GRBmsg(env, message)
    ccall((:GRBmsg, libgurobi), Cvoid, (Ptr{GRBenv}, Ptr{Cchar}), env, message)
end

function GRBgetlogfile(env, logfileP)
    ccall((:GRBgetlogfile, libgurobi), Cint, (Ptr{GRBenv}, Ptr{Ptr{Cint}}), env, logfileP)
end

function GRBsetlogfile(env, logfile)
    ccall((:GRBsetlogfile, libgurobi), Cint, (Ptr{GRBenv}, Ptr{Cint}), env, logfile)
end

function GRBgetintparam(env, paramname, valueP)
    ccall((:GRBgetintparam, libgurobi), Cint, (Ptr{GRBenv}, Ptr{Cchar}, Ptr{Cint}), env, paramname, valueP)
end

function GRBgetdblparam(env, paramname, valueP)
    ccall((:GRBgetdblparam, libgurobi), Cint, (Ptr{GRBenv}, Ptr{Cchar}, Ptr{Cdouble}), env, paramname, valueP)
end

function GRBgetstrparam(env, paramname, valueP)
    ccall((:GRBgetstrparam, libgurobi), Cint, (Ptr{GRBenv}, Ptr{Cchar}, Ptr{Cchar}), env, paramname, valueP)
end

function GRBgetlongstrparam(env, paramname, valueP, size, requiredlenP)
    ccall((:GRBgetlongstrparam, libgurobi), Cint, (Ptr{GRBenv}, Ptr{Cchar}, Ptr{Cchar}, Cint, Ptr{Cint}), env, paramname, valueP, size, requiredlenP)
end

function GRBgetintparaminfo(env, paramname, valueP, minP, maxP, defP)
    ccall((:GRBgetintparaminfo, libgurobi), Cint, (Ptr{GRBenv}, Ptr{Cchar}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}), env, paramname, valueP, minP, maxP, defP)
end

function GRBgetdblparaminfo(env, paramname, valueP, minP, maxP, defP)
    ccall((:GRBgetdblparaminfo, libgurobi), Cint, (Ptr{GRBenv}, Ptr{Cchar}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}), env, paramname, valueP, minP, maxP, defP)
end

function GRBgetstrparaminfo(env, paramname, valueP, defP)
    ccall((:GRBgetstrparaminfo, libgurobi), Cint, (Ptr{GRBenv}, Ptr{Cchar}, Ptr{Cchar}, Ptr{Cchar}), env, paramname, valueP, defP)
end

function GRBgetparamflags(env, parname, valueP)
    ccall((:GRBgetparamflags, libgurobi), Cint, (Ptr{GRBenv}, Ptr{Cchar}, Ptr{UInt32}), env, parname, valueP)
end

function GRBsetparam(env, paramname, value)
    ccall((:GRBsetparam, libgurobi), Cint, (Ptr{GRBenv}, Ptr{Cchar}, Ptr{Cchar}), env, paramname, value)
end

function GRBsetintparam(env, paramname, value)
    ccall((:GRBsetintparam, libgurobi), Cint, (Ptr{GRBenv}, Ptr{Cchar}, Cint), env, paramname, value)
end

function GRBsetdblparam(env, paramname, value)
    ccall((:GRBsetdblparam, libgurobi), Cint, (Ptr{GRBenv}, Ptr{Cchar}, Cdouble), env, paramname, value)
end

function GRBsetstrparam(env, paramname, value)
    ccall((:GRBsetstrparam, libgurobi), Cint, (Ptr{GRBenv}, Ptr{Cchar}, Ptr{Cchar}), env, paramname, value)
end

function GRBgetparamtype(env, paramname)
    ccall((:GRBgetparamtype, libgurobi), Cint, (Ptr{GRBenv}, Ptr{Cchar}), env, paramname)
end

function GRBresetparams(env)
    ccall((:GRBresetparams, libgurobi), Cint, (Ptr{GRBenv},), env)
end

function GRBcopyparams(dest, src)
    ccall((:GRBcopyparams, libgurobi), Cint, (Ptr{GRBenv}, Ptr{GRBenv}), dest, src)
end

function GRBwriteparams(env, filename)
    ccall((:GRBwriteparams, libgurobi), Cint, (Ptr{GRBenv}, Ptr{Cchar}), env, filename)
end

function GRBreadparams(env, filename)
    ccall((:GRBreadparams, libgurobi), Cint, (Ptr{GRBenv}, Ptr{Cchar}), env, filename)
end

function GRBgetnumparams(env)
    ccall((:GRBgetnumparams, libgurobi), Cint, (Ptr{GRBenv},), env)
end

function GRBgetparamname(env, i, paramnameP)
    ccall((:GRBgetparamname, libgurobi), Cint, (Ptr{GRBenv}, Cint, Ptr{Ptr{Cchar}}), env, i, paramnameP)
end

function GRBgetnumattributes(model)
    ccall((:GRBgetnumattributes, libgurobi), Cint, (Ptr{GRBmodel},), model)
end

function GRBgetattrname(model, i, attrnameP)
    ccall((:GRBgetattrname, libgurobi), Cint, (Ptr{GRBmodel}, Cint, Ptr{Ptr{Cchar}}), model, i, attrnameP)
end

function GRBloadenv(envP, logfilename)
    ccall((:GRBloadenv, libgurobi), Cint, (Ptr{Ptr{GRBenv}}, Ptr{Cchar}), envP, logfilename)
end

function GRBemptyenv(envP)
    ccall((:GRBemptyenv, libgurobi), Cint, (Ptr{Ptr{GRBenv}},), envP)
end

function GRBstartenv(env)
    ccall((:GRBstartenv, libgurobi), Cint, (Ptr{GRBenv},), env)
end

function GRBloadenvadv(envP, logfilename, apitype, major, minor, tech, server, router, password, group, priority, idletimeout, accessid, secretkey, cb, usrdata, logcb)
    ccall((:GRBloadenvadv, libgurobi), Cint, (Ptr{Ptr{GRBenv}}, Ptr{Cchar}, Cint, Cint, Cint, Cint, Ptr{Cchar}, Ptr{Cchar}, Ptr{Cchar}, Ptr{Cchar}, Cint, Cint, Ptr{Cchar}, Ptr{Cchar}, Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}), envP, logfilename, apitype, major, minor, tech, server, router, password, group, priority, idletimeout, accessid, secretkey, cb, usrdata, logcb)
end

function GRBloadclientenv(envP, logfilename, computeserver, router, password, group, CStlsinsecure, priority, timeout)
    ccall((:GRBloadclientenv, libgurobi), Cint, (Ptr{Ptr{GRBenv}}, Ptr{Cchar}, Ptr{Cchar}, Ptr{Cchar}, Ptr{Cchar}, Ptr{Cchar}, Cint, Cint, Cdouble), envP, logfilename, computeserver, router, password, group, CStlsinsecure, priority, timeout)
end

function GRBloadclientenvadv(envP, logfilename, computeserver, router, password, group, CStlsinsecure, priority, timeout, apitype, major, minor, tech, cb, usrdata)
    ccall((:GRBloadclientenvadv, libgurobi), Cint, (Ptr{Ptr{GRBenv}}, Ptr{Cchar}, Ptr{Cchar}, Ptr{Cchar}, Ptr{Cchar}, Ptr{Cchar}, Cint, Cint, Cdouble, Cint, Cint, Cint, Cint, Ptr{Cvoid}, Ptr{Cvoid}), envP, logfilename, computeserver, router, password, group, CStlsinsecure, priority, timeout, apitype, major, minor, tech, cb, usrdata)
end

function GRBloadcloudenv(envP, logfilename, accessID, secretKey, pool, priority)
    ccall((:GRBloadcloudenv, libgurobi), Cint, (Ptr{Ptr{GRBenv}}, Ptr{Cchar}, Ptr{Cchar}, Ptr{Cchar}, Ptr{Cchar}, Cint), envP, logfilename, accessID, secretKey, pool, priority)
end

function GRBloadcloudenvadv(envP, logfilename, accessID, secretKey, pool, priority, apitype, major, minor, tech, cb, usrdata)
    ccall((:GRBloadcloudenvadv, libgurobi), Cint, (Ptr{Ptr{GRBenv}}, Ptr{Cchar}, Ptr{Cchar}, Ptr{Cchar}, Ptr{Cchar}, Cint, Cint, Cint, Cint, Cint, Ptr{Cvoid}, Ptr{Cvoid}), envP, logfilename, accessID, secretKey, pool, priority, apitype, major, minor, tech, cb, usrdata)
end

function GRBgetenv(model)
    ccall((:GRBgetenv, libgurobi), Ptr{GRBenv}, (Ptr{GRBmodel},), model)
end

function GRBgetconcurrentenv(model, num)
    ccall((:GRBgetconcurrentenv, libgurobi), Ptr{GRBenv}, (Ptr{GRBmodel}, Cint), model, num)
end

function GRBdiscardconcurrentenvs(model)
    ccall((:GRBdiscardconcurrentenvs, libgurobi), Cvoid, (Ptr{GRBmodel},), model)
end

function GRBgetmultiobjenv(model, num)
    ccall((:GRBgetmultiobjenv, libgurobi), Ptr{GRBenv}, (Ptr{GRBmodel}, Cint), model, num)
end

function GRBdiscardmultiobjenvs(model)
    ccall((:GRBdiscardmultiobjenvs, libgurobi), Cvoid, (Ptr{GRBmodel},), model)
end

function GRBgettuneenv(model, num)
    ccall((:GRBgettuneenv, libgurobi), Ptr{GRBenv}, (Ptr{GRBmodel}, Cint), model, num)
end

function GRBdiscardtuneenvs(model)
    ccall((:GRBdiscardtuneenvs, libgurobi), Cvoid, (Ptr{GRBmodel},), model)
end

function GRBreleaselicense(env)
    ccall((:GRBreleaselicense, libgurobi), Cvoid, (Ptr{GRBenv},), env)
end

function GRBfreeenv(env)
    ccall((:GRBfreeenv, libgurobi), Cvoid, (Ptr{GRBenv},), env)
end

function GRBgeterrormsg(env)
    ccall((:GRBgeterrormsg, libgurobi), Ptr{Cchar}, (Ptr{GRBenv},), env)
end

function GRBgetmerrormsg(model)
    ccall((:GRBgetmerrormsg, libgurobi), Ptr{Cchar}, (Ptr{GRBmodel},), model)
end

function GRBgetcommstats(env, recvtimeP, recvbytesP, recvmsgsP, sendtimeP, sendbytesP, sendmsgsP)
    ccall((:GRBgetcommstats, libgurobi), Cvoid, (Ptr{GRBenv}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}), env, recvtimeP, recvbytesP, recvmsgsP, sendtimeP, sendbytesP, sendmsgsP)
end

function GRBversion(majorP, minorP, technicalP)
    ccall((:GRBversion, libgurobi), Cvoid, (Ptr{Cint}, Ptr{Cint}, Ptr{Cint}), majorP, minorP, technicalP)
end

function GRBplatform()
    ccall((:GRBplatform, libgurobi), Ptr{Cchar}, ())
end

function GRBlisttokens()
    ccall((:GRBlisttokens, libgurobi), Cint, ())
end

function GRBprinttuneparams()
    ccall((:GRBprinttuneparams, libgurobi), Cvoid, ())
end

function GRBtunemodel(model)
    ccall((:GRBtunemodel, libgurobi), Cint, (Ptr{GRBmodel},), model)
end

function GRBtunemodels(nummodels, models, ignore)
    ccall((:GRBtunemodels, libgurobi), Cint, (Cint, Ptr{Ptr{GRBmodel}}, Ptr{GRBmodel}), nummodels, models, ignore)
end

function GRBgettuneresult(model, i)
    ccall((:GRBgettuneresult, libgurobi), Cint, (Ptr{GRBmodel}, Cint), model, i)
end

function GRBgettunelog(model, i, logP)
    ccall((:GRBgettunelog, libgurobi), Cint, (Ptr{GRBmodel}, Cint, Ptr{Ptr{Cchar}}), model, i, logP)
end

function GRBtunemodeladv(model, ignore)
    ccall((:GRBtunemodeladv, libgurobi), Cint, (Ptr{GRBmodel}, Ptr{GRBmodel}), model, ignore)
end

function GRBsortIDi(len, ind, val)
    ccall((:GRBsortIDi, libgurobi), Cvoid, (Cint, Ptr{Cint}, Ptr{Cdouble}), len, ind, val)
end

function GRBabortbatch(batch)
    ccall((:GRBabortbatch, libgurobi), Cint, (Ptr{GRBbatch},), batch)
end

function GRBdiscardbatch(batch)
    ccall((:GRBdiscardbatch, libgurobi), Cint, (Ptr{GRBbatch},), batch)
end

function GRBretrybatch(batch)
    ccall((:GRBretrybatch, libgurobi), Cint, (Ptr{GRBbatch},), batch)
end

function GRBfreebatch(batch)
    ccall((:GRBfreebatch, libgurobi), Cint, (Ptr{GRBbatch},), batch)
end

function GRBgetbatch(env, batchID, batchP)
    ccall((:GRBgetbatch, libgurobi), Cint, (Ptr{GRBenv}, Ptr{Cchar}, Ptr{Ptr{GRBbatch}}), env, batchID, batchP)
end

function GRBgetbatchjsonsolution(batch, jsonsolP)
    ccall((:GRBgetbatchjsonsolution, libgurobi), Cint, (Ptr{GRBbatch}, Ptr{Ptr{Cchar}}), batch, jsonsolP)
end

function GRBgetbatchintattr(batch, attrname, valueP)
    ccall((:GRBgetbatchintattr, libgurobi), Cint, (Ptr{GRBbatch}, Ptr{Cchar}, Ptr{Cint}), batch, attrname, valueP)
end

function GRBgetbatchstrattr(batch, attrname, valueP)
    ccall((:GRBgetbatchstrattr, libgurobi), Cint, (Ptr{GRBbatch}, Ptr{Cchar}, Ptr{Ptr{Cchar}}), batch, attrname, valueP)
end

function GRBgetbatchattrname(env, n, attrnameP)
    ccall((:GRBgetbatchattrname, libgurobi), Cint, (Ptr{GRBenv}, Cint, Ptr{Ptr{Cchar}}), env, n, attrnameP)
end

function GRBgetbatchattrflags(batch, attrname, flagsP)
    ccall((:GRBgetbatchattrflags, libgurobi), Cint, (Ptr{GRBbatch}, Ptr{Cchar}, Ptr{UInt32}), batch, attrname, flagsP)
end

function GRBgetbatchattrinfo(batch, attrname, datatypeP, settableP)
    ccall((:GRBgetbatchattrinfo, libgurobi), Cint, (Ptr{GRBbatch}, Ptr{Cchar}, Ptr{Cint}, Ptr{Cint}), batch, attrname, datatypeP, settableP)
end

function GRBupdatebatch(batch)
    ccall((:GRBupdatebatch, libgurobi), Cint, (Ptr{GRBbatch},), batch)
end

function GRBwritebatchjsonsolution(batch, filename)
    ccall((:GRBwritebatchjsonsolution, libgurobi), Cint, (Ptr{GRBbatch}, Ptr{Cchar}), batch, filename)
end

function GRBgetnumbatchattributes(env)
    ccall((:GRBgetnumbatchattributes, libgurobi), Cint, (Ptr{GRBenv},), env)
end

function GRBgetbatchenv(batch)
    ccall((:GRBgetbatchenv, libgurobi), Ptr{GRBenv}, (Ptr{GRBbatch},), batch)
end

function GRBfree(ptr)
    ccall((:GRBfree, libgurobi), Cvoid, (Ptr{Cvoid},), ptr)
end

function GRBsync(model)
    ccall((:GRBsync, libgurobi), Cint, (Ptr{GRBmodel},), model)
end

function GRBpingserver(server, password)
    ccall((:GRBpingserver, libgurobi), Cint, (Ptr{Cchar}, Ptr{Cchar}), server, password)
end

function GRBprefetchattr(model, attrname)
    ccall((:GRBprefetchattr, libgurobi), Cint, (Ptr{GRBmodel}, Ptr{Cchar}), model, attrname)
end
