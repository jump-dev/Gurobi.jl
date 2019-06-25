# Variables in Gurobi model

# variable kinds

const GRB_CONTINUOUS = convert(Cchar, 'C')
const GRB_BINARY     = convert(Cchar, 'B')
const GRB_INTEGER    = convert(Cchar, 'I')

# add_var!
function add_var!(model::Model, numnz::Integer, vind::Vector, vval::Vector{Float64}, c::Float64, lb::Float64, ub::Float64, vtype::Cchar)
    if checkvalue(c, GRB_INFINITY)
        _objwarning(c)
    end
    if checkvalue(lb, GRB_BOUNDMAX) || checkvalue(ub, GRB_BOUNDMAX)
        _boundwarning(lb, ub)
    end
    ret = @grb_ccall(addvar, Cint, (
        Ptr{Cvoid},    # model
        Cint,         # numnz
        Ptr{Cint},    # vind
        Ptr{Float64}, # vval
        Float64,      # obj
        Float64,      # lb
        Float64,      # ub
        UInt8,        # vtype
        Ptr{UInt8}    # name
        ),
        model, numnz, ivec(vind.-1), vval, c, lb, ub, vtype, C_NULL)

    if ret != 0
        throw(GurobiError(model.env, ret))
    end
    nothing
end

function add_var!(model::Model, vtype::Cchar, c::Float64, lb::Float64, ub::Float64)
    if checkvalue(c, GRB_INFINITY)
        _objwarning(c)
    end
    if checkvalue(lb, GRB_BOUNDMAX) || checkvalue(ub, GRB_BOUNDMAX)
        _boundwarning(lb, ub)
    end
    ret = @grb_ccall(addvar, Cint, (
        Ptr{Cvoid},    # model
        Cint,         # numnz
        Ptr{Cint},    # vind
        Ptr{Float64}, # vval
        Float64,      # obj
        Float64,      # lb
        Float64,      # ub
        UInt8,        # vtype
        Ptr{UInt8}    # name
        ),
        model, 0, C_NULL, C_NULL, c, lb, ub, vtype, C_NULL)

    if ret != 0
        throw(GurobiError(model.env, ret))
    end
    nothing
end

add_var!(model::Model, vtype::GChars, c::Real, lb::Real, ub::Real) = add_var!(model, cchar(vtype), Float64(c), Float64(lb), Float64(ub))
add_var!(model::Model, vtype::GChars, c::Real) = add_var!(model, vtype, c, -Inf, Inf)

add_cvar!(model::Model, c::Real, lb::Real, ub::Real) = add_var!(model, GRB_CONTINUOUS, c, lb, ub)
add_cvar!(model::Model, c::Real) = add_cvar!(model, c, -Inf, Inf)

add_bvar!(model::Model, c::Real) = add_var!(model, GRB_BINARY, c, 0., 1.)

add_ivar!(model::Model, c::Real, lb::Real, ub::Real) = add_var!(model, GRB_INTEGER, c, lb, ub)
add_ivar!(model::Model, c::Real) = add_ivar!(model, c, -Inf, Inf)


# add_vars!

function add_vars!(model::Model, vtypes::CVec, c::FVec, lb::FVec, ub::FVec)
    if checkvalue(c, GRB_INFINITY)
        _objwarning(c)
    end
    if checkvalue(lb, GRB_BOUNDMAX) || checkvalue(ub, GRB_BOUNDMAX)
        _boundwarning(lb, ub)
    end
    # check dimensions
    n = length(vtypes)
    _chklen(c, n)
    _chklen(lb, n)
    _chklen(ub, n)

    # main call
    ret = @grb_ccall(addvars, Cint, (
        Ptr{Cvoid},  # model
        Cint,       # numvars
        Cint,       # numnz
        Ptr{Cint},  # vbeg
        Ptr{Cint},  # vind
        Ptr{Float64}, # vval
        Ptr{Float64}, # obj
        Ptr{Float64}, # lb
        Ptr{Float64}, # ub
        Ptr{Cchar},   # vtypes
        Ptr{Ptr{UInt8}}, # varnames
        ),
        model, n, 0, C_NULL, C_NULL, C_NULL, c, lb, ub, vtypes, C_NULL)

    if ret != 0
        throw(GurobiError(model.env, ret))
    end
    nothing
end

function add_vars!(model::Model, vtypes::GCharOrVec, c::Vector, lb::Bounds, ub::Bounds)
    n = length(c)
    add_vars!(model, cvecx(vtypes, n), fvec(c), fvecx(lb, n), fvecx(ub, n))
end

add_cvars!(model::Model, c::Vector, lb::Bounds, ub::Bounds) = add_vars!(model, GRB_CONTINUOUS, c, lb, ub)
add_cvars!(model::Model, c::Vector) = add_cvars!(model, c, -Inf, Inf)

add_bvars!(model::Model, c::Vector) = add_vars!(model, GRB_BINARY, c, 0, 1)

add_ivars!(model::Model, c::Vector, lb::Bounds, ub::Bounds) = add_vars!(model, GRB_INTEGER, c, lb, ub)
add_ivars!(model::Model, c::Vector) = add_ivars!(model, GRB_INTEGER, c, -Inf, Inf)

del_vars!(model::Model, idx::T) where {T<:Real} = del_vars!(model, Cint[idx])
del_vars!(model::Model, idx::Vector{T}) where {T<:Real} = del_vars!(model, convert(Vector{Cint},idx))
function del_vars!(model::Model, idx::Vector{Cint})
    numdel = length(idx)
    ret = @grb_ccall(delvars, Cint, (
                     Ptr{Cvoid},
                     Cint,
                     Ptr{Cint}),
                     model, convert(Cint,numdel), idx.-Cint(1))
    if ret != 0
        throw(GurobiError(model.env, ret))
    end
end

function get_vars(model::Model, start::Integer, len::Integer)
# http://www.gurobi.com/documentation/8.1/refman/c_grbgetvars.html

    # NR_CONTRAINTS
    #--------------
    m = num_constrs(model)
    # NR_VARIABLES
    #-------------
    n = num_vars(model)

    # INPUT VALIDATION
    #-----------------
    @assert start > 0               "Indexing in Julia starts from 1."
    @assert start <= n              string("Index out of bounds: There are only ", n, " variables attached to this model.")
    @assert len > 0                 "At least one variable must be selected; len > 0."
    @assert start * len <= m * n    "Maximal amount of possible non-zero elements surpassed."

    # FUNCTION CALLS
    #---------------
    numnzP = Ref{Cint}()
    vbeg = Array{Cint}(undef, len)

    ret = @grb_ccall(getvars, Cint, (
                        Ptr{Cvoid},     #Model
                        Ptr{Cint},      #numnzP
                        Ptr{Cint},      #vbeg
                        Ptr{Cint},      #vind
                        Ptr{Cdouble},   #vval
                        Cint,           #start
                        Cint),          #len
                        model, numnzP, vbeg, C_NULL, C_NULL, Cint(start - 1), Cint(len))
    if ret != 0
        throw(GurobiError(model.env, ret))
    end

    nnz = numnzP[]
    vind = Array{Cint}(undef, nnz)
    vval = Array{Cdouble}(undef, nnz)

    ret = @grb_ccall(getvars, Cint, (
                        Ptr{Cvoid},     #Model
                        Ptr{Cint},      #numnzP
                        Ptr{Cint},      #vbeg
                        Ptr{Cint},      #vind
                        Ptr{Cdouble},   #vval
                        Cint,           #start
                        Cint),          #len
                        model, numnzP, vbeg, vind, vval, Cint(start - 1), Cint(len))

    if ret != 0
        throw(GurobiError(model.env, ret))
    end

    # ADJUSTING INDICES TO JULIA'S INDEXING.
    #---------------------------------------
    for i in 1:size(vbeg, 1)
        vbeg[i] += 1
    end

    for i in 1:size(vind, 1)
        vind[i] += 1
    end

    # SPARSE ARRAY
    #-------------
    push!(vbeg, nnz)
    I = Array{Int64}(undef, nnz)
    J = Array{Int64}(undef, nnz)
    V = Array{Float64}(undef, nnz)
    for i in 1:length(vbeg) - 1
        for j in vbeg[i]:vbeg[i + 1]
            I[j] = vind[j]
            J[j] = i + start - 1
            V[j] = vval[j]
        end
    end

    pop!(vbeg)

    # return vbeg, vind, vval
    return SparseArrays.sparse(I, J, V, m, n)
end
