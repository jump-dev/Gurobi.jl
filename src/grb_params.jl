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
    "DegenMoves",
    "Disconnected",
    "MinRelNodes",
    "MIPFocus",
    "NodeMethod",
    "NoRelHeuristic",
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
    "ProjImpliedCuts",
    "MIPSepCuts",
    "MIRCuts",
    "StrongCGCuts",
    "ModKCuts",
    "ZeroHalfCuts",
    "NetworkCuts",
    "SubMIPCuts",
    "InfProofCuts",
    "CutAggPasses",
    "CutPasses",
    "GomoryPasses",
    "WorkerPort",
    "Aggregate",
    "AggFill",
    "ConcurrentMIP",
    "ConcurrentMIPJobs", # deprecated after Gurobi7
    "ConcurrentJobs",
    "DisplayInterval",
    "DistributedMIPJobs",
    "DualReductions",
    "IISMethod",
    "InfUnbdInfo",
    "LazyConstraints",
    "LogToConsole",
    "MIQCPMethod",
    "NumericFocus",
    "NonBlocking", # deprecated after Gurobi7
    "OutputFlag",
    "PreCrush",
    "PreDepRow",
    "PreDual",
    "PrePasses",
    "PreQLinearize",
    "Presolve",
    "PreSparsify",
    "PreMIQCPForm",
    "QCPDual",
    "Record",
    "Seed",
    "Threads",
    "TuneResults",
    "TuneCriterion",
    "TuneTrials",
    "TuneOutput",
    "TuneJobs",
    "PreMIQPMethod", # deprecated after Gurobi7
    "UpdateMode",
    "ObjNumber",
    "MultiObjMethod",
    "MultiObjPre",
    "PoolSolutions",
    "PoolSearchMode"]

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
    "BestObjStop",
    "BestBdStop",
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
    "TuneTimeLimit",
    "PoolGap"]

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

function get_int_param(env::Env, name::String)
    @assert isascii(name)
    a = Ref{Cint}()
    ret = @grb_ccall(getintparam, Cint, (Ptr{Cvoid}, Ptr{Cchar}, Ptr{Cint}),
        env, name, a)
    if ret != 0
        throw(GurobiError(env, ret))
    end
    convert(Int, a[])
end

function get_dbl_param(env::Env, name::String)
    @assert isascii(name)
    a = Ref{Float64}()
    ret = @grb_ccall(getdblparam, Cint, (Ptr{Cvoid}, Ptr{Cchar}, Ptr{Float64}),
        env, name, a)
    if ret != 0
        throw(GurobiError(env, ret))
    end
    a[]::Float64
end

function get_str_param(env::Env, name::String)
    @assert isascii(name)
    buf = Array{Cchar}(undef, GRB_MAX_STRLEN)
    ret = @grb_ccall(getstrparam, Cint, (Ptr{Cvoid}, Ptr{Cchar}, Ptr{Cchar}),
        env, name, buf)
    if ret != 0
        throw(GurobiError(env, ret))
    end
    unsafe_string(pointer(buf))
end


function set_int_param!(env::Env, name::String, v::Integer)
    @assert isascii(name)
    ret = @grb_ccall(setintparam, Cint, (Ptr{Cvoid}, Ptr{Cchar}, Cint),
        env, name, v)
    if ret != 0
        throw(GurobiError(env, ret))
    end
end

function set_dbl_param!(env::Env, name::String, v::Real)
    @assert isascii(name)
    ret = @grb_ccall(setdblparam, Cint, (Ptr{Cvoid}, Ptr{Cchar}, Float64),
        env, name, v)
    if ret != 0
        throw(GurobiError(env, ret))
    end
end

function set_str_param!(env::Env, name::String, v::String)
    @assert isascii(name)
    @assert isascii(v)
    ret = @grb_ccall(setstrparam, Cint, (Ptr{Cvoid}, Ptr{Cchar}, Ptr{Cchar}),
        env, name, v)
    if ret != 0
        throw(GurobiError(env, ret))
    end
end

# for existing models

function get_int_param(m::Model, name::String)

    @assert isascii(name)
    modenv = @grb_ccall(getenv, Ptr{Cvoid}, (Ptr{Cvoid},), m)
    @assert modenv != C_NULL
    a = Ref{Cint}()
    ret = @grb_ccall(getintparam, Cint, (Ptr{Cvoid}, Ptr{Cchar}, Ptr{Cint}),
        modenv, name, a)
    if ret != 0
        throw(GurobiError(Env(modenv), ret))
    end
    convert(Int, a[])
end

function get_dbl_param(m::Model, name::String)

    @assert isascii(name)
    modenv = @grb_ccall(getenv, Ptr{Cvoid}, (Ptr{Cvoid},), m)
    @assert modenv != C_NULL
    a = Ref{Float64}()
    ret = @grb_ccall(getdblparam, Cint, (Ptr{Cvoid}, Ptr{Cchar}, Ptr{Float64}),
        modenv, name, a)
    if ret != 0
        throw(GurobiError(Env(modenv), ret))
    end
    a[]::Float64
end

function get_str_param(m::Model, name::String)

    @assert isascii(name)
    modenv = @grb_ccall(getenv, Ptr{Cvoid}, (Ptr{Cvoid},), m)
    @assert modenv != C_NULL
    buf = Array{Cchar}(undef, GRB_MAX_STRLEN)
    ret = @grb_ccall(getstrparam, Cint, (Ptr{Cvoid}, Ptr{Cchar}, Ptr{Cchar}),
        modenv, name, buf)
    if ret != 0
        throw(GurobiError(Env(modenv), ret))
    end
    unsafe_string(pointer(buf))
end


function set_int_param!(m::Model, name::String, v::Integer)

    @assert isascii(name)
    modenv = @grb_ccall(getenv, Ptr{Cvoid}, (Ptr{Cvoid},), m)
    @assert modenv != C_NULL
    ret = @grb_ccall(setintparam, Cint, (Ptr{Cvoid}, Ptr{Cchar}, Cint),
        modenv, name, v)
    if ret != 0
        throw(GurobiError(Env(modenv), ret))
    end
end

function set_dbl_param!(m::Model, name::String, v::Real)

    @assert isascii(name)
    modenv = @grb_ccall(getenv, Ptr{Cvoid}, (Ptr{Cvoid},), m)
    @assert modenv != C_NULL
    ret = @grb_ccall(setdblparam, Cint, (Ptr{Cvoid}, Ptr{Cchar}, Float64),
        modenv, name, v)
    if ret != 0
        throw(GurobiError(Env(modenv), ret))
    end
end

function set_str_param!(m::Model, name::String, v::String)

    @assert isascii(name)
    @assert isascii(v)
    modenv = @grb_ccall(getenv, Ptr{Cvoid}, (Ptr{Cvoid},), m)
    @assert modenv != C_NULL
    ret = @grb_ccall(setstrparam, Cint, (Ptr{Cvoid}, Ptr{Cchar}, Ptr{Cchar}),
        modenv, name, v)
    if ret != 0
        throw(GurobiError(Env(modenv), ret))
    end
end

# higher level functions

function getparam(env::Union{Env,Model}, name::String)
    @assert isascii(name)
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

function setparam!(env::Union{Env,Model}, name::String, v)
    @assert isascii(name)
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
