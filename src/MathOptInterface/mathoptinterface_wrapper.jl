export GurobiOptimizer

# using Gurobi
# const GRB = Gurobi
using MathOptInterface
const MOI = MathOptInterface
using LinQuadOptInterface
const LQOI = LinQuadOptInterface

const SUPPORTED_OBJECTIVES = [
    LQOI.Linear,
    LQOI.Quad
]
const SUPPORTED_CONSTRAINTS = [
    (LQOI.Linear, LQOI.EQ),
    (LQOI.Linear, LQOI.LE),
    (LQOI.Linear, LQOI.GE),
    # (Linear, IV),
    (LQOI.Quad, LQOI.EQ),
    (LQOI.Quad, LQOI.LE),
    (LQOI.Quad, LQOI.GE),
    (LQOI.SinVar, LQOI.EQ),
    (LQOI.SinVar, LQOI.LE),
    (LQOI.SinVar, LQOI.GE),
    (LQOI.SinVar, LQOI.IV),
    (LQOI.SinVar, MOI.ZeroOne),
    (LQOI.SinVar, MOI.Integer),
    (LQOI.VecVar, LQOI.SOS1),
    (LQOI.VecVar, LQOI.SOS2),
    (LQOI.VecVar, MOI.Nonnegatives),
    (LQOI.VecVar, MOI.Nonpositives),
    (LQOI.VecVar, MOI.Zeros),
    (LQOI.VecLin, MOI.Nonnegatives),
    (LQOI.VecLin, MOI.Nonpositives),
    (LQOI.VecLin, MOI.Zeros)
]

mutable struct GurobiOptimizer <: LQOI.LinQuadOptimizer
    LQOI.@LinQuadOptimizerBase
    env
    params::Dict{String,Any}
    GurobiOptimizer(::Void) = new()
end

LQOI.LinQuadModel(::Type{GurobiOptimizer},env) = Model(env::Env,"defaultname")

function GurobiOptimizer(;kwargs...)

    env = Env()
    m = GurobiOptimizer(nothing)
    m.env = env
    m.params = Dict{String,Any}()
    MOI.empty!(m)
    for (name,value) in kwargs
        m.params[string(name)] = value
        setparam!(m.inner, string(name), value)
    end
    return m
end

function MOI.empty!(m::GurobiOptimizer)
    MOI.empty!(m,m.env)
    for (name,value) in m.params
        setparam!(m.inner, name, value)
    end
end

LQOI.lqs_supported_constraints(s::GurobiOptimizer) = SUPPORTED_CONSTRAINTS
LQOI.lqs_supported_objectives(s::GurobiOptimizer) = SUPPORTED_OBJECTIVES
#=
    inner wrapper
=#

#=
    Main
=#

# LinQuadSolver # Abstract type
# done above

# LQOI.lqs_setparam!(env, name, val)
# TODO fix this one
LQOI.lqs_setparam!(m::GurobiOptimizer, name, val) = setparam!(m.inner, string(name), val)

# LQOI.lqs_setlogfile!(env, path)
# TODO fix this one
LQOI.lqs_setlogfile!(m::GurobiOptimizer, path) = setlogfile(m.env, path::String)

# LQOI.lqs_getprobtype(m)
# TODO - consider removing, apparently useless

#=
    Constraints
=#

cintvec(v::Vector) = convert(Vector{Int32}, v)

_getsense(m::GurobiOptimizer, ::MOI.EqualTo{Float64}) = Cchar('=')
_getsense(m::GurobiOptimizer, ::MOI.LessThan{Float64}) = Cchar('<')
_getsense(m::GurobiOptimizer, ::MOI.GreaterThan{Float64}) = Cchar('>')
_getsense(m::GurobiOptimizer, ::MOI.Zeros)        = Cchar('=')
_getsense(m::GurobiOptimizer, ::MOI.Nonpositives) = Cchar('<')
_getsense(m::GurobiOptimizer, ::MOI.Nonnegatives) = Cchar('>')
_getboundsense(m::GurobiOptimizer, ::MOI.Nonpositives) = Cchar('>')
_getboundsense(m::GurobiOptimizer, ::MOI.Nonnegatives) = Cchar('<')

# LQOI.lqs_chgbds!(m, colvec, valvec, sensevec)
# TODO - improve single type
function LQOI.lqs_chgbds!(instance::GurobiOptimizer, colvec, valvec, sensevec)
    lb_len = count(x->x==Cchar('L'), sensevec)
    LB_val = Array{Float64}(0)
    sizehint!(LB_val, lb_len)
    LB_col = Array{Cint}(0)
    sizehint!(LB_col, lb_len)

    ub_len = count(x->x==Cchar('U'), sensevec)
    UB_val = Array{Float64}(0)
    sizehint!(UB_val, ub_len)
    UB_col = Array{Cint}(0)
    sizehint!(UB_col, ub_len)

    for i in eachindex(valvec)
        if sensevec[i] == Cchar('L')
            push!(LB_col, colvec[i])
            push!(LB_val, valvec[i])
        elseif sensevec[i] == Cchar('U')
            push!(UB_col, colvec[i])
            push!(UB_val, valvec[i])
        end
    end

    if lb_len > 0
        set_dblattrlist!(instance.inner, "LB", LB_col, LB_val)
    end

    if ub_len > 0
        set_dblattrlist!(instance.inner, "UB", UB_col, UB_val)
    end

    update_model!(instance.inner)
    nothing
end


# LQOI.lqs_getlb(m, col)
LQOI.lqs_getlb(instance::GurobiOptimizer, col) = (update_model!(instance.inner);get_dblattrlist( instance.inner, "LB", ivec(col))[1])
# LQOI.lqs_getub(m, col)
LQOI.lqs_getub(instance::GurobiOptimizer, col) = get_dblattrlist( instance.inner, "UB", ivec(col))[1]

# LQOI.lqs_getnumrows(m)
LQOI.lqs_getnumrows(instance::GurobiOptimizer) = num_constrs(instance.inner)

# LQOI.lqs_addrows!(m, rowvec, colvec, coefvec, sensevec, rhsvec)
LQOI.lqs_addrows!(instance::GurobiOptimizer, rowvec, colvec, coefvec, sensevec, rhsvec) = (add_constrs!(instance.inner, rowvec, colvec, coefvec, sensevec, rhsvec);update_model!(instance.inner))

# LQOI.lqs_getrhs(m, rowvec)
LQOI.lqs_getrhs(instance::GurobiOptimizer, row) = get_dblattrlist( instance.inner, "RHS", ivec(row))[1]

# colvec, coef = LQOI.lqs_getrows(m, rowvec)
# TODO improve
function LQOI.lqs_getrows(instance::GurobiOptimizer, idx)
    A = get_constrs(instance.inner, idx, 1)'
    return A.rowval-1, A.nzval
end

# LQOI.lqs_getcoef(m, row, col) #??
# TODO improve
function LQOI.lqs_getcoef(instance::GurobiOptimizer, row, col) #??
    return getcoeff(model::Model, row::Integer, col::Integer)
    # A = get_rows(m, row, row)'
    # cols = A.rowval
    # vals = A.nzval

    # pos = findfirst(cols, col)
    # if pos > 0
    #     return vals[pos]
    # else
    #     return 0.0
    # end
end

# LQOI.lqs_chgcoef!(m, row, col, coef)
# TODO SPLIT THIS ONE
function LQOI.lqs_chgcoef!(instance::GurobiOptimizer, row, col, coef)
    if row == 0
        set_dblattrlist!(instance.inner, "Obj", Cint[col], Float64[coef])
    elseif col == 0
        set_dblattrlist!(instance.inner, "RHS", Cint[row], Float64[coef])
    else
        chg_coeffs!(instance.inner, row, col, coef)
        #TODO fix this function in gurobi
    end
end

# LQOI.lqs_delrows!(m, row, row)
LQOI.lqs_delrows!(instance::GurobiOptimizer, rowbeg, rowend) = del_constrs!(instance.inner, cintvec(collect(rowbeg:rowend)))

# LQOI.lqs_chgctype!(m, colvec, typevec)
# TODO fix types
LQOI.lqs_chgctype!(instance::GurobiOptimizer, colvec, typevec) = set_charattrlist!(instance.inner, "VType", ivec(colvec), cvec(typevec))

# LQOI.lqs_chgsense!(m, rowvec, sensevec)
# TODO fix types
LQOI.lqs_chgsense!(instance::GurobiOptimizer, rowvec, sensevec) = set_charattrlist!(instance.inner, "Sense", ivec(rowvec), cvec(sensevec))

const VAR_TYPE_MAP = Dict{Symbol,Cchar}(
    :CONTINUOUS => Cchar('C'),
    :INTEGER => Cchar('I'),
    :BINARY => Cchar('B')
)
LQOI.lqs_vartype_map(m::GurobiOptimizer) = VAR_TYPE_MAP

# LQOI.lqs_addsos(m, colvec, valvec, typ)
LQOI.lqs_addsos!(instance::GurobiOptimizer, colvec, valvec, typ) = (add_sos!(instance.inner, typ, colvec, valvec);update_model!(instance.inner))
# LQOI.lqs_delsos(m, idx, idx)
LQOI.lqs_delsos!(instance::GurobiOptimizer, idx1, idx2) = (del_sos!(instance.inner, cintvec(collect(idx1:idx2)));update_model!(instance.inner))

const SOS_TYPE_MAP = Dict{Symbol,Symbol}(
    :SOS1 => :SOS1,#Cchar('1'),
    :SOS2 => :SOS2#Cchar('2')
)
LQOI.lqs_sertype_map(m::GurobiOptimizer) = SOS_TYPE_MAP

# LQOI.lqs_getsos(m, idx)
# TODO improve getting processes
function LQOI.lqs_getsos(instance::GurobiOptimizer, idx)
    A, types = get_sos_matrix(instance.inner)
    line = A[idx,:] #sparse vec
    cols = line.nzind
    vals = line.nzval
    typ = types[idx] == Cint(1) ? :SOS1 : :SOS2
    return cols, vals, typ
end

# LQOI.lqs_getnumqconstrs(m)
LQOI.lqs_getnumqconstrs(instance::GurobiOptimizer) = num_qconstrs(instance.inner)

# LQOI.lqs_addqconstr(m, cols,coefs,rhs,sense, I,J,V)
#   NOTE:
# LQOI assumes 0.5 x' Q x, but Gurobi requires x' Q x so we multiply V by 0.5
LQOI.lqs_addqconstr!(instance::GurobiOptimizer, cols,coefs,rhs,sense, I,J,V) = add_qconstr!(instance.inner, cols, coefs, I, J, 0.5 * V, sense, rhs)

# LQOI.lqs_chgrngval
LQOI.lqs_chgrngval!(instance::GurobiOptimizer, rows, vals) = chg_rhsrange!(instance.inner, cintvec(rows), -vals)

const CTR_TYPE_MAP = Dict{Symbol,Cchar}(
    :RANGE => Cchar('R'),
    :LOWER => Cchar('L'),
    :UPPER => Cchar('U'),
    :EQUALITY => Cchar('E')
)
LQOI.lqs_ctrtype_map(m::GurobiOptimizer) = CTR_TYPE_MAP

#=
    Objective
=#

# LQOI.lqs_copyquad(m, intvec,intvec, floatvec) #?
function LQOI.lqs_copyquad!(instance::GurobiOptimizer, I, J, V)
    delq!(instance.inner)
    for i in eachindex(V)
        if I[i] == J[i]
            V[i] /= 2
        end
    end
    add_qpterms!(instance.inner, I, J, V)
    return nothing
end

# LQOI.lqs_chgobj(m, colvec,coefvec)
function LQOI.lqs_chgobj!(instance::GurobiOptimizer, colvec, coefvec)
    nvars = num_vars(instance.inner)
    obj = zeros(Float64, nvars)

    for i in eachindex(colvec)
        obj[colvec[i]] = coefvec[i]
    end

    set_dblattrarray!(instance.inner, "Obj", 1, num_vars(instance.inner), obj)
    update_model!(instance.inner)
    nothing
end

# LQOI.lqs_chgobjsen(m, symbol)
# TODO improve min max names
function LQOI.lqs_chgobjsen!(instance::GurobiOptimizer, symbol)
    if symbol in [:minimize, :Min]
        set_sense!(instance.inner, :minimize)
    else
        set_sense!(instance.inner, :maximize)
    end
    update_model!(instance.inner)
end


# LQOI.lqs_getobj(instance.inner)
LQOI.lqs_getobj(instance::GurobiOptimizer) = get_dblattrarray( instance.inner, "Obj", 1, num_vars(instance.inner)   )

# lqs_getobjsen(m)
function LQOI.lqs_getobjsen(instance::GurobiOptimizer)
    s = model_sense(instance.inner)
    if s in [:maximize, :Max]
        return MOI.MaxSense
    else
        return MOI.MinSense
    end
end

#=
    Variables
=#

# LQOI.lqs_getnumcols(m)
LQOI.lqs_getnumcols(instance::GurobiOptimizer) = (update_model!(instance.inner); num_vars(instance.inner))

# LQOI.lqs_newcols!(m, int)
LQOI.lqs_newcols!(instance::GurobiOptimizer, int) = (add_cvars!(instance.inner, zeros(int));update_model!(instance.inner))

# LQOI.lqs_delcols!(m, col, col)
LQOI.lqs_delcols!(instance::GurobiOptimizer, col, col2) = (del_vars!(instance.inner, col);update_model!(instance.inner))

# LQOI.lqs_addmipstarts(m, colvec, valvec)
function LQOI.lqs_addmipstarts!(instance::GurobiOptimizer, colvec, valvec)
    x = zeros(num_vars(instance.inner))
    for i in eachindex(colvec)
        x[colvec[i]] = valvec[i]
    end
    loadbasis(instance.inner, x)
end
#=
    Solve
=#

# LQOI.lqs_mipopt!(m)
LQOI.lqs_mipopt!(instance::GurobiOptimizer) = LQOI.lqs_lpopt!(instance)

# LQOI.lqs_qpopt!(m)
LQOI.lqs_qpopt!(instance::GurobiOptimizer) = LQOI.lqs_lpopt!(instance)

# LQOI.lqs_lpopt!(m)
LQOI.lqs_lpopt!(instance::GurobiOptimizer) = (update_model!(instance.inner);optimize(instance.inner))

# LQOI.lqs_terminationstatus(m)
function LQOI.lqs_terminationstatus(instance::GurobiOptimizer)

    stat = get_status(instance.inner)

    if stat == :loaded
        return MOI.OtherError
    elseif stat == :optimal
        return MOI.Success
    elseif stat == :infeasible
        if hasdualray(instance)
            return MOI.Success
        else
            return MOI.InfeasibleNoResult
        end
    elseif stat == :inf_or_unbd
        return MOI.InfeasibleOrUnbounded
    elseif stat == :unbounded
        if hasprimalray(instance)
            return MOI.Success
        else
            return MOI.UnboundedNoResult
        end
    elseif stat == :cutoff
        return MOI.ObjectiveLimit
    elseif stat == :iteration_limit
        return MOI.IterationLimit
    elseif stat == :node_limit
        return MOI.NodeLimit
    elseif stat == :time_limit
        return MOI.TimeLimit
    elseif stat == :solution_limit
        return MOI.SolutionLimit
    elseif stat == :interrupted
        return MOI.Interrupted
    elseif stat == :numeric
        return MOI.NumericalError
    elseif stat == :suboptimal
        return MOI.OtherLimit
    elseif stat == :inprogress
        return MOI.OtherError
    elseif stat == :user_obj_limit
        return MOI.ObjectiveLimit
    end
    return MOI.OtherError
end


function LQOI.lqs_primalstatus(instance::GurobiOptimizer)

    stat = get_status(instance.inner)

    if stat == :optimal
        return MOI.FeasiblePoint
    elseif stat == :solution_limit
        return MOI.FeasiblePoint
    elseif stat in [:inf_or_unbd, :unbounded] && hasprimalray(instance)
        return MOI.InfeasibilityCertificate
    elseif stat == :suboptimal
        return MOI.FeasiblePoint
    else
        return MOI.UnknownResultStatus
    end
end
function LQOI.lqs_dualstatus(instance::GurobiOptimizer)
    stat = get_status(instance.inner)

    if is_mip(instance.inner) || is_qcp(instance.inner)
        return MOI.UnknownResultStatus
    else
        if stat == :optimal
            return MOI.FeasiblePoint
        elseif stat == :solution_limit
            return MOI.FeasiblePoint
        elseif stat in [:inf_or_unbd, :infeasible] && hasdualray(instance)
            return MOI.InfeasibilityCertificate
        elseif stat == :suboptimal
            return MOI.FeasiblePoint
        else
            return MOI.UnknownResultStatus
        end
    end
end


# LQOI.lqs_getx!(m, place)
LQOI.lqs_getx!(instance::GurobiOptimizer, place) = get_dblattrarray!(place, instance.inner, "X", 1)

# LQOI.lqs_getax!(m, place)
function LQOI.lqs_getax!(instance::GurobiOptimizer, place)
    get_dblattrarray!(place, instance.inner, "Slack", 1)
    rhs = get_dblattrarray(instance.inner, "RHS", 1, num_constrs(instance.inner))
    for i in eachindex(place)
        place[i] = -place[i]+rhs[i]
    end
    nothing
end
# LQOI.lqs_getdj!(m, place)
LQOI.lqs_getdj!(instance::GurobiOptimizer, place) = get_dblattrarray!(place, instance.inner, "RC", 1)

# LQOI.lqs_getpi!(m, place)
LQOI.lqs_getpi!(instance::GurobiOptimizer, place) = get_dblattrarray!(place, instance.inner, "Pi", 1)

# LQOI.lqs_getobjval(m)
LQOI.lqs_getobjval(instance::GurobiOptimizer) = get_objval(instance.inner)

# LQOI.lqs_getbestobjval(m)
LQOI.lqs_getbestobjval(instance::GurobiOptimizer) = get_objval(instance.inner)

# LQOI.lqs_getmiprelgap(m)
function LQOI.lqs_getmiprelgap(instance::GurobiOptimizer)
    L = get_objval(instance.inner)
    U = get_objbound(instance.inner)
    return abs(U-L)/U
end

# LQOI.lqs_getitcnt(m)
LQOI.lqs_getitcnt(instance::GurobiOptimizer)  = get_iter_count(instance.inner)

# LQOI.lqs_getbaritcnt(m)
LQOI.lqs_getbaritcnt(instance::GurobiOptimizer) = get_barrier_iter_count(instance.inner)

# LQOI.lqs_getnodecnt(m)
LQOI.lqs_getnodecnt(instance::GurobiOptimizer) = get_node_count(instance.inner)

# LQOI.lqs_dualfarkas(m, place)
LQOI.lqs_dualfarkas!(instance::GurobiOptimizer, place) = get_dblattrarray!(place, instance.inner, "FarkasDual", 1)

function hasdualray(instance::GurobiOptimizer)
    try
        get_dblattrarray(instance.inner, "FarkasDual", 1, num_constrs(instance.inner))
        return true
    catch
        return false
    end
end

# LQOI.lqs_getray(m, place)
LQOI.lqs_getray!(instance::GurobiOptimizer, place) = get_dblattrarray!(place, instance.inner, "UnbdRay", 1)

function hasprimalray(instance::GurobiOptimizer)
    try
        get_dblattrarray(instance.inner, "UnbdRay", 1, num_vars(instance.inner))
        return true
    catch
        return false
    end
end

MOI.free!(m::GurobiOptimizer) = free_model(m.inner)

"""
    writeproblem(m: :MOI.AbstractOptimizer, filename::String)
Writes the current problem data to the given file.
Supported file types are solver-dependent.
"""
writeproblem(m::GurobiOptimizer, filename::String, flags::String="") = write_model(m.inner, filename)


# blocked
MOI.addconstraint!(m::GurobiOptimizer, func::LQOI.Linear, set::LQOI.IV) = error("not supported")
MOI.addconstraints!(m::GurobiOptimizer, func::Vector{LQOI.Linear}, set::Vector{LQOI.IV}) = error("not supported")

MOI.canget(m::GurobiOptimizer, any, c::LQOI.LCI{LQOI.IV}) = false
MOI.canmodifyconstraint(m::GurobiOptimizer, c::LQOI.LCI{LQOI.IV}, chg) = false
MOI.candelete(m::GurobiOptimizer, c::LQOI.LCI{LQOI.IV}) = false
