# Variables in Gurobi model

# variable kinds

const GRB_CONTINUOUS = convert(Cchar, 'C')
const GRB_BINARY     = convert(Cchar, 'B')
const GRB_INTEGER    = convert(Cchar, 'I')

# add_var!

function add_var!(model::Model, vtype::Cchar, c::Float64, lb::Float64, ub::Float64)
    
    ret = @grb_ccall(addvar, Cint, (
        Ptr{Void},    # model
        Cint,         # numnz
        Ptr{Cint},    # vind
        Ptr{Float64}, # vval
        Float64,      # obj
        Float64,      # lb
        Float64,      # ub
        Uint8,        # vtype
        Ptr{Uint8}    # name
        ), 
        model, 0, C_NULL, C_NULL, c, lb, ub, vtype, C_NULL)
        
    if ret != 0
        throw(GurobiError(model.env, ret))
    end
    nothing
end

add_var!(model::Model, vtype::Cchar, c::Float64) = add_var!(model, vtype, c, -Inf, Inf)

add_cvar!(model::Model, c::Float64, lb::Float64, ub::Float64) = add_var!(model, GRB_CONTINUOUS, c, lb, ub)
add_cvar!(model::Model, c::Float64) = add_cvar!(model, c, -Inf, Inf)

add_bvar!(model::Model, c::Float64) = add_var!(model, GRB_BINARY, c, 0., 1.)

add_ivar!(model::Model, c::Float64, lb::Real, ub::Real) = add_var!(model, GRB_INTEGER, c, float64(lb), float64(ub))
add_ivar!(model::Model, c::Float64) = add_ivar!(model, c, -Inf, Inf)


# add_vars!

function add_vars!(model::Model, vtypes::Union(Cchar, Vector{Cchar}), c::Vector{Float64}, lb::Bounds, ub::Bounds)
    
    # check dimensions
    
    n = length(c)    
    dims_ok::Bool = true
    
    if isa(vtypes, Vector) && length(vtypes) != n
        dims_ok = false
    end

    if isa(lb, Vector) && length(lb) != n
        dims_ok = false
    end
    
    if isa(ub, Vector) && length(lb) != n
        dims_ok = false
    end
    
    if !dims_ok
        throw(ArgumentError("Inconsistent dimensions."))
    end

    # prepare vtypes
    
    vtypes_ptr::Ptr{Cchar} = C_NULL
    
    if isa(vtypes, Cchar)
        if vtypes != GRB_CONTINUOUS
            vtypes = fill(vtypes, n)
            vtypes_ptr = pointer(vtypes)
        end
    else
        vtypes_ptr = pointer(vtypes)
    end

    # prepare bounds
    
    lb_ptr::Ptr{Float64} = C_NULL
    ub_ptr::Ptr{Float64} = C_NULL
    
    if lb == nothing
        lb = fill(-Inf, n)
        lb_ptr = pointer(lb)
    elseif isa(lb, Float64)
        if lb != 0.
            lb = fill(lb, n)
            lb_ptr = pointer(lb)
        end
    else
        lb_ptr = pointer(lb)
    end
    
    if isa(ub, Float64)
        if !(ub == Inf)
            ub = fill(ub, n)
            ub_ptr = pointer(ub)
        end
    elseif isa(ub, Vector)
        ub_ptr = pointer(ub)
    end

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
        Ptr{Ptr{Uint8}}, # varnames
        ), 
        model, n, 0, C_NULL, C_NULL, C_NULL, 
        c, lb_ptr, ub_ptr, vtypes_ptr, C_NULL)
        
    if ret != 0
        throw(GurobiError(model.env, ret))
    end
    nothing
end

add_cvars!(model::Model, c::Vector{Float64}, lb::Bounds, ub::Bounds) = add_vars!(model, GRB_CONTINUOUS, c, lb, ub)
add_cvars!(model::Model, c::Vector{Float64}) = add_cvars!(model, c, nothing, nothing)

add_bvars!(model::Model, c::Vector{Float64}) = add_vars!(model, GRB_BINARY, c, 0., 1.)

add_ivars!(model::Model, c::Vector{Float64}, lb::Bounds, ub::Bounds) = add_vars!(model, GRB_INTEGER, c, lb, ub)
add_ivars!(model::Model, c::Vector{Float64}) = add_ivars!(model, GRB_INTEGER, c, nothing, nothing) 

