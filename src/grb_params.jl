# Gurobi solver parameters

# This list was obtained through AWK (with Gurobi 5.6)
# command:
#
#    grep GRB_INT_PAR gurobi_c.h | awk '{ printf("%s,\n", $3) }'
#

const GRB_INT_PARAMS =[ 
    "BarIterLimit",
    "SolutionLimit",
    "Method",
    "ScaleFlag",
    "SimplexPricing",
    "Quad",
    "NormAdjust",
    "Sifting",
    "SiftMethod",
    "BarCorrectors",
    "BarHomogeneous",
    "BarOrder",
    "Crossover",
    "CrossoverBasis",
    "BranchDir",
    "Disconnected",
    "MinRelNodes",
    "MIPFocus",
    "NodeMethod",
    "PumpPasses",
    "RINS",
    "SubMIPNodes",
    "Symmetry",
    "VarBranch",
    "SolutionNumber",
    "ZeroObjNodes",
    "Cuts",
    "CliqueCuts",
    "CoverCuts",
    "FlowCoverCuts",
    "FlowPathCuts",
    "GUBCoverCuts",
    "ImpliedCuts",
    "MIPSepCuts",
    "MIRCuts",
    "ModKCuts",
    "ZeroHalfCuts",
    "NetworkCuts",
    "SubMIPCuts",
    "CutAggPasses",
    "CutPasses",
    "GomoryPasses",
    "Aggregate",
    "AggFill",
    "ConcurrentMIP",
    "ConcurrentMIPJobs",
    "DisplayInterval",
    "DualReductions",
    "IISMethod",
    "InfUnbdInfo",
    "LazyConstraints",
    "LogToConsole",
    "MIQCPMethod",
    "NumericFocus",
    "NonBlocking",
    "OutputFlag",
    "PreCrush",
    "PreDepRow",
    "PreDual",
    "PrePasses",
    "PreQLinearize",
    "Presolve",
    "PreSparsify",
    "QCPDual",
    "Seed",
    "Threads",
    "TuneResults",
    "TuneTrials",
    "TuneOutput",
    "TuneJobs",
    "PreMIQPMethod"]

# This list was obtained through AWK (with Gurobi 5.6)
# command:
#
#    grep GRB_DBL_PAR gurobi_c.h | awk '{ printf("%s,\n", $3) }'
#
const GRB_DBL_PARAMS = [
    "Cutoff",
    "IterationLimit",
    "NodeLimit",
    "TimeLimit",
    "FeasibilityTol",
    "IntFeasTol",
    "MarkowitzTol",
    "MIPGap",
    "MIPGapAbs",
    "OptimalityTol",
    "PSDTol",
    "PerturbValue",
    "ObjScale",
    "BarConvTol",
    "BarQCPConvTol",
    "Heuristics",
    "ImproveStartGap",
    "ImproveStartTime",
    "ImproveStartNodes",
    "NodefileStart",
    "FeasRelaxBigM",
    "PreSOS1BigM",
    "PreSOS2BigM",
    "TuneTimeLimit"]

# This list was obtained through AWK (with Gurobi 5.6)
# command:
#
#    grep GRB_STR_PAR gurobi_c.h | awk '{ printf("%s,\n", $3) }'
#
const GRB_STR_PARAMS = [
    "NodefileDir",
    "ServerPool",
    "ServerPassword",
    "LogFile",
    "ResultFile",
    "Dummy"]


#################################################
#
#  Parameter getter & setter
#
#################################################

const GRB_MAX_STRLEN = 512   # gurobi.c define this value as 512

# lower-level functions

function get_int_param(env::Env, name::ASCIIString)
    a = Array(Cint, 1)
    ret = @grb_ccall(getintparam, Cint, (Ptr{Void}, Ptr{Cchar}, Ptr{Cint}), 
        env, name, a)
    if ret != 0
        throw(GurobiError(env, ret))
    end
    convert(Int, a[1])
end

function get_dbl_param(env::Env, name::ASCIIString)
    a = Array(Float64, 1)
    ret = @grb_ccall(getdblparam, Cint, (Ptr{Void}, Ptr{Cchar}, Ptr{Float64}), 
        env, name, a)
    if ret != 0
        throw(GurobiError(env, ret))
    end
    a[1]::Float64
end

function get_str_param(env::Env, name::ASCIIString)
    buf = Array(Cchar, GRB_MAX_STRLEN)
    ret = @grb_ccall(getstrparam, Cint, (Ptr{Void}, Ptr{Cchar}, Ptr{Cchar}), 
        env, name, buf)
    if ret != 0
        throw(GurobiError(env, ret))
    end
    bytestring(pointer(buf))
end


function set_int_param!(env::Env, name::ASCIIString, v::Integer)
    ret = @grb_ccall(setintparam, Cint, (Ptr{Void}, Ptr{Cchar}, Cint), 
        env, name, v)
    if ret != 0
        throw(GurobiError(env, ret))
    end
end

function set_dbl_param!(env::Env, name::ASCIIString, v::Real)
    ret = @grb_ccall(setdblparam, Cint, (Ptr{Void}, Ptr{Cchar}, Float64), 
        env, name, v)
    if ret != 0
        throw(GurobiError(env, ret))
    end
end

function set_str_param!(env::Env, name::ASCIIString, v::ASCIIString)
    ret = @grb_ccall(setstrparam, Cint, (Ptr{Void}, Ptr{Cchar}, Ptr{Cchar}), 
        env, name, v)
    if ret != 0
        throw(GurobiError(env, ret))        
    end
end

# for existing models

function get_int_param(m::Model, name::ASCIIString)

    modenv = @grb_ccall(getenv, Ptr{Void}, (Ptr{Void},), m)
    @assert modenv != C_NULL
    a = Array(Cint, 1)
    ret = @grb_ccall(getintparam, Cint, (Ptr{Void}, Ptr{Cchar}, Ptr{Cint}),
        modenv, name, a)
    if ret != 0
        throw(GurobiError(Env(modenv), ret))
    end
    convert(Int, a[1])
end

function get_dbl_param(m::Model, name::ASCIIString)

    modenv = @grb_ccall(getenv, Ptr{Void}, (Ptr{Void},), m)
    @assert modenv != C_NULL
    a = Array(Float64, 1)
    ret = @grb_ccall(getdblparam, Cint, (Ptr{Void}, Ptr{Cchar}, Ptr{Float64}),
        modenv, name, a)
    if ret != 0
        throw(GurobiError(Env(modenv), ret))
    end
    a[1]::Float64
end

function get_str_param(m::Model, name::ASCIIString)

    modenv = @grb_ccall(getenv, Ptr{Void}, (Ptr{Void},), m)
    @assert modenv != C_NULL
    buf = Array(Cchar, GRB_MAX_STRLEN)
    ret = @grb_ccall(getstrparam, Cint, (Ptr{Void}, Ptr{Cchar}, Ptr{Cchar}),
        modenv, name, buf)
    if ret != 0
        throw(GurobiError(Env(modenv), ret))
    end
    bytestring(pointer(buf))
end


function set_int_param!(m::Model, name::ASCIIString, v::Integer)

    modenv = @grb_ccall(getenv, Ptr{Void}, (Ptr{Void},), m)
    @assert modenv != C_NULL
    ret = @grb_ccall(setintparam, Cint, (Ptr{Void}, Ptr{Cchar}, Cint),
        modenv, name, v)
    if ret != 0
        throw(GurobiError(Env(modenv), ret))
    end
end

function set_dbl_param!(m::Model, name::ASCIIString, v::Real)

    modenv = @grb_ccall(getenv, Ptr{Void}, (Ptr{Void},), m)
    @assert modenv != C_NULL
    ret = @grb_ccall(setdblparam, Cint, (Ptr{Void}, Ptr{Cchar}, Float64),
        modenv, name, v)
    if ret != 0
        throw(GurobiError(Env(modenv), ret))
    end
end

function set_str_param!(m::Model, name::ASCIIString, v::ASCIIString)

    modenv = @grb_ccall(getenv, Ptr{Void}, (Ptr{Void},), m)
    @assert modenv != C_NULL
    ret = @grb_ccall(setstrparam, Cint, (Ptr{Void}, Ptr{Cchar}, Ptr{Cchar}),
        modenv, name, v)
    if ret != 0
        throw(GurobiError(Env(modenv), ret))
    end
end

# higher level functions

function getparam(env::Union{Env,Model}, name::ASCIIString)
    if name in GRB_INT_PARAMS
        return get_int_param(env, name)
    elseif name in GRB_DBL_PARAMS
        return get_dbl_param(env, name)
    elseif name in GRB_STR_PARAMS
        return get_str_param(env, name)
    else
        error("Unrecognized parameter name: $(name).")
    end
end

function setparam!(env::Union{Env,Model}, name::ASCIIString, v)
    if name in GRB_INT_PARAMS
        set_int_param!(env, name, v)
    elseif name in GRB_DBL_PARAMS
        set_dbl_param!(env, name, v)
    elseif name in GRB_STR_PARAMS
        set_str_param!(env, name, v)
    else
        error("Unrecognized parameter name: $(name).")
    end
end

function setparams!(env::Union{Env,Model}; args...)
    for (name, v) in args
        setparam!(env, string(name), v)
    end
end

