# Variables in Gurobi model

# variable kinds

const GRB_CONTINUOUS = convert(Cchar, 'C')
const GRB_BINARY     = convert(Cchar, 'B')
const GRB_INTEGER    = convert(Cchar, 'I')

# add_var!
function add_var!(model::Model, numnz::Integer, vind::Vector, vval::Vector{Float64}, obj::Float64, lb::Float64, ub::Float64, vtype::Cchar)
    ret = @grb_ccall(addvar, Cint, (
        Ptr{Void},    # model
        Cint,         # numnz
        Ptr{Cint},    # vind
        Ptr{Float64}, # vval
        Float64,      # obj
        Float64,      # lb
        Float64,      # ub
        UInt8,        # vtype
        Ptr{UInt8}    # name
        ), 
        model, numnz, ivec(vind.-1), vval, obj, lb, ub, vtype, C_NULL)
        
    if ret != 0
        throw(GurobiError(model.env, ret))
    end
    nothing
end

function add_var!(model::Model, vtype::Cchar, c::Float64, lb::Float64, ub::Float64)    
    ret = @grb_ccall(addvar, Cint, (
        Ptr{Void},    # model
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

function add_var!(model::Model, vtype::GChars, c::Real, lb::Real, ub::Real)
    add_var!(model, cchar(vtype), Float64(c), Float64(lb), Float64(ub))
end
add_var!(model::Model, vtype::GChars, c::Real) = add_var!(model, vtype, c, -Inf, Inf)

add_cvar!(model::Model, c::Real, lb::Real, ub::Real) = add_var!(model, GRB_CONTINUOUS, c, lb, ub)
add_cvar!(model::Model, c::Real) = add_cvar!(model, c, -Inf, Inf)

add_bvar!(model::Model, c::Real) = add_var!(model, GRB_BINARY, c, 0., 1.)

add_ivar!(model::Model, c::Real, lb::Real, ub::Real) = add_var!(model, GRB_INTEGER, c, lb, ub)
add_ivar!(model::Model, c::Real) = add_ivar!(model, c, -Inf, Inf)


# add_vars!

function add_vars!(model::Model, vtypes::CVec, c::FVec, lb::FVec, ub::FVec)
    
    # check dimensions
    n = length(vtypes)    
    _chklen(c, n)
    _chklen(lb, n)
    _chklen(ub, n)

    # main call
    ret = @grb_ccall(addvars, Cint, (
        Ptr{Void},  # model
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

