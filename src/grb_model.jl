# Gurobi model

type Model
    env::Env
    ptr_model::Ptr{Void}
    
    function Model(env::Env, p::Ptr{Void})
        model = new(env, p)
        finalizer(model, free_model)
        model
    end
end

convert(ty::Type{Ptr{Void}}, model::Model) = model.ptr_model::Ptr{Void}

function free_model(model::Model)
    if model.ptr_model != C_NULL
        @grb_ccall(freemodel, Void, (Ptr{Void},), model.ptr_model)
        model.ptr_model = C_NULL
    end
end

function copy(model::Model)
    pm::Ptr{Void} = C_NULL
    if model.ptr_model != C_NULL
        pm = @grb_ccall(copymodel, Ptr{Void}, (Ptr{Void},), model.ptr_model)
        if pm == C_NULL
            error("Failed to copy a Gurobi model.")
        end
    end
    Model(model.env, pm)
end

function update_model!(model::Model)
    @assert model.ptr_model != C_NULL
    ret = @grb_ccall(updatemodel, Cint, (Ptr{Void},), model.ptr_model)
    if ret != 0
        throw(GurobiError(model.env, ret))
    end
    nothing
end

function reset_model!(model::Model)
    @assert model.ptr_model != C_NULL
    ret = @grb_ccall(resetmodel, Cint, (Ptr{Void},), model.ptr_model)
    if ret != 0
        throw(GurobiError(model.env, ret))
    end
    nothing
end

# read / write file

function read_model(env::Env, filename::ASCIIString)
    a = Array(Ptr{Void}, 1)
    ret = @grb_ccall(readmodel, Cint, 
        (Ptr{Void}, Ptr{Uint8}, Ptr{Ptr{Void}}), 
        env, filename, a)
    if ret != 0
        throw(GurobiError(env, ret))
    end
    Model(env, a[1])
end

function write_model(model::Model, filename::ASCIIString)
    ret = @grb_ccall(write, Cint, (Ptr{Void}, Ptr{Uint8}), 
        model.ptr_model, filename)
    if ret != 0
        throw(GurobiError(model.env, ret))
    end
    nothing
end



#################################################
#
#  model attributes
#
#################################################

function get_int_attr(model::Model, name::ASCIIString)
    a = Array(Cint, 1)
    ret = @grb_ccall(getintattr, Cint, 
        (Ptr{Void}, Ptr{Uint8}, Ptr{Cint}),
        model, name, a);
    if ret != 0
        throw(GurobiError(model.env, ret))
    end
    convert(Int, a[1])
end

function get_dbl_attr(model::Model, name::ASCIIString)
    a = Array(Float64, 1)
    ret = @grb_ccall(getdblattr, Cint, 
        (Ptr{Void}, Ptr{Uint8}, Ptr{Float64}),
        model, name, a);
    if ret != 0
        throw(GurobiError(model.env, ret))
    end
    a[1]::Float64
end

function get_str_attr(model::Model, name::ASCIIString)
    a = Array(Ptr{Uint8}, 1)
    ret = @grb_ccall(getstrattr, Cint, 
        (Ptr{Void}, Ptr{Uint8}, Ptr{Ptr{Uint8}}), 
        model, name, a)
    if ret != 0
        throw(GurobiError(model.env, ret))
    end
    bytestring(a[1])
end

function get_dbl_attrarray(model::Model, name::ASCIIString, start::Integer, len::Integer)
    # start is one-based
    
    r = Array(Float64, len)
    ret = @grb_ccall(getdblattrarray, Cint, 
        (Ptr{Void}, Ptr{Uint8}, Cint, Cint, Ptr{Float64}), 
        model, name, start - 1, len, pointer(r))
    if ret != 0
        throw(GurobiError(model.env, ret))
    end
    r
end

function set_int_attr!(model::Model, name::ASCIIString, v::Integer)
    ret = @grb_ccall(setintattr, Cint, 
        (Ptr{Void}, Ptr{Uint8}, Cint), model, name, v)
    if ret != 0
        throw(GurobiError(model.env, ret))
    end
    nothing
end

function set_dbl_attr!(model::Model, name::ASCIIString, v::Real)
    ret = @grb_ccall(setdblattr, Cint, 
        (Ptr{Void}, Ptr{Uint8}, Float64), model, name, v)
    if ret != 0
        throw(GurobiError(model.env, ret))
    end
    nothing
end

function set_str_attr!(model::Model, name::ASCIIString, v::ASCIIString)
    ret = @grb_ccall(setstrattr, Cint, 
        (Ptr{Void}, Ptr{Uint8}, Ptr{Uint8}), model, name, v)
    if ret != 0
        throw(GurobiError(model.env, ret))
    end
    nothing
end

function set_char_attr_array!(model::Model, name::ASCIIString, start::Integer, len::Integer, values::Vector{Char})
    values = convert(Vector{Cchar},values)
    ret = @grb_ccall(setcharattrarray, Cint, 
        (Ptr{Void}, Ptr{Uint8}, Cint, Cint, Ptr{Cchar}), model, name, start-1, len, values)
    if ret != 0
        throw(GurobiError(model.env, ret))
    end
    nothing
end

macro grb_int_attr(fun, attrname)
    @eval $(fun)(model::Model) = get_int_attr(model, $(attrname))
end

macro grb_dbl_attr(fun, attrname)
    @eval $(fun)(model::Model) = get_dbl_attr(model, $(attrname))
end

macro grb_str_attr(fun, attrname)
    @eval $(fun)(model::Model) = get_str_attr(model, $(attrname))
end

@grb_int_attr num_vars     "NumVars"
@grb_int_attr num_constrs  "NumConstrs"
@grb_int_attr num_sos      "NumSOS"
@grb_int_attr num_qconstrs "NumQConstrs"
@grb_int_attr num_cnzs     "NumNZs"
@grb_int_attr num_qnzs     "NumQNZs"
@grb_int_attr num_qcnzs    "NumCQNZs"

@grb_str_attr model_name  "ModelName"

name(model::Model) = model_name(model)
sense(model::Model) = get_int_attr(model, "ModelSense") > 0 ? (:minimize) : (:maximize)

is_mip(model::Model) = get_int_attr(model, "IsMIP") != 0
is_qp(model::Model)  = get_int_attr(model, "IsQP") != 0
is_qcp(model::Model) = get_int_attr(model, "IsQCP") != 0

function model_type(model::Model) 
    is_qp(model)  ? (:QP)  :
    is_qcp(model) ? (:QCP) : (:LP)
end

function set_sense!(model::Model, sense::Symbol)
    v = sense == :maximize ? -1 :
        sense == :minimize ? 1 : 
        throw(ArgumentError("Invalid model sense."))
    
    set_int_attr!(model, "ModelSense", v)
end


# Show

function show(io::IO, model::Model)
    if model.ptr_model != C_NULL
        println(io, "Gurobi Model: $(name(model))")
        if is_mip(model)
            println(io, "    type   : $(model_type(model)) (MIP)")
        else
            println(io, "    type   : $(model_type(model))")
        end
        println(io, "    sense  : $(sense(model))")
        println(io, "    number of variables             = $(num_vars(model))")
        println(io, "    number of linear constraints    = $(num_constrs(model))")
        println(io, "    number of quadratic constraints = $(num_qconstrs(model))")
        println(io, "    number of sos constraints       = $(num_sos(model))")
        println(io, "    number of non-zero coeffs       = $(num_cnzs(model))")
        println(io, "    number of non-zero qp-terms     = $(num_qnzs(model))")
    else
        println(io, "Gurobi Model: NULL")
    end
end


#################################################
#
#  model construction
#
#################################################

typealias Bounds Union(Nothing, Float64, Vector{Float64})

function gurobi_model(env::Env, name::ASCIIString)
    @assert is_valid(env)
    
    a = Array(Ptr{Void}, 1)
    ret = @grb_ccall(newmodel, Cint, (
        Ptr{Void},  # env
        Ptr{Ptr{Void}},  # modelp
        Ptr{Uint8},  # name
        Cint,  # numvars
        Ptr{Float64}, # obj coeffs
        Ptr{Float64}, # lbounds
        Ptr{Float64}, # ubounds
        Ptr{Uint8},   # var types,
        Ptr{Ptr{Uint8}} # varnames
        ), 
        env, a, name, 0, 
        C_NULL, C_NULL, C_NULL, C_NULL, C_NULL
    )
    
    if ret != 0
        throw(GurobiError(env, ret))
    end
    
    Model(env, a[1])
end


function gurobi_model(env::Env, name::ASCIIString, sense::Symbol)
    model = gurobi_model(env, name)
    if sense != :minimize
        set_sense!(model, sense)
    end
    model
end


# add variables

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

# add_constr

add_constr!(model::Model, inds::Vector, coeffs::Vector{Float64}, rel::Char, rhs::Float64) =
    add_constr!(model, convert(Vector{Cint},inds), coeffs, rel, rhs)

function add_constr!(model::Model, inds::Vector{Cint}, coeffs::Vector{Float64}, rel::Char, rhs::Float64)
    inds = inds - 1
    if !isempty(inds)
        ret = @grb_ccall(addconstr, Cint, (
            Ptr{Void},    # model
            Cint,         # numnz
            Ptr{Cint},    # cind
            Ptr{Float64}, # cvals
            Cchar,        # sense
            Float64,      # rhs
            Ptr{Uint8}    # name
            ), 
            model, length(inds), inds, coeffs, rel, rhs, C_NULL)
        if ret != 0
            throw(GurobiError(model.env, ret))
        end
    end
    nothing 
end

function add_constr!(model::Model, coeffs::Vector{Float64}, rel::Char, rhs::Float64)
    inds = convert(Vector{Cint}, find(coeffs)) 
    vals = coeffs[inds]
    add_constr!(model, inds, vals, rel, rhs)
end

# add_constrs!

function add_constrs!(model::Model, cbegins::Vector{Cint}, inds::Vector{Cint}, coeffs::Vector{Float64}, 
    senses::Vector{Cchar}, rhs::Vector{Float64})
        
    m = length(cbegins)
    nnz = length(inds)
    
    if !(m == length(senses) == length(rhs) && nnz == length(coeffs))
        throw(ArgumentError("Incompatible dimensions."))
    end 
        
    if m > 0 && nnz > 0
        ret = @grb_ccall(addconstrs, Cint, (
            Ptr{Void},    # model
            Cint,         # num constraints
            Cint,         # num non-zeros
            Ptr{Cint},    # cbeg
            Ptr{Cint},    # cind
            Ptr{Float64}, # cval
            Ptr{Cchar},   # sense
            Ptr{Float64}, # rhs
            Ptr{Uint8}    # names
            ), 
            model, m, nnz, cbegins - 1, inds - 1, coeffs, 
            senses, rhs, C_NULL)
        
        if ret != 0
            throw(GurobiError(model.env, ret))
        end
    end
    nothing
end

function add_constrs!(
    model::Model, cbegins::Vector{Cint}, inds::Vector{Cint}, coeffs::Vector{Float64}, 
    rel::Char, rhs::Vector{Float64})
    add_constrs!(model, cbegins, inds, coeffs, fill(convert(Cchar, rel), length(cbegins)), rhs)
end


# add_qterms!

add_qpterms!(model::Model, qr::Vector, qc::Vector, qv::Vector{Float64}) =
    add_qpterms!(model, convert(Vector{Cint}, qr), convert(Vector{Cint}, qc),qv)

function add_qpterms!(model::Model, qr::Vector{Cint}, qc::Vector{Cint}, qv::Vector{Float64})
    nnz = length(qr)
    if !(nnz == length(qc) == length(qv))
        throw(ArgumentError("Inconsistent dimensions."))
    end
    
    if nnz > 0
        ret = @grb_ccall(addqpterms, Cint, (
            Ptr{Void},    # model
            Cint,         # nnz
            Ptr{Cint},    # qrow
            Ptr{Cint},    # qcol
            Ptr{Float64}, # qval
            ), 
            model, nnz, qr-1, qc-1, qv)
            
        if ret != 0
            throw(GurobiError(model.env, ret))
        end 
    end
    nothing
end

# add_qconstr!

add_qconstr!(model::Model, lind::Vector, lval::Vector, qr::Vector, qc::Vector,
    qv::Vector{Float64}, rel::Char, rhs::Float64) =
    add_qconstr!(model, convert(Vector{Cint},lind), convert(Vector{Float64}, lval),
        convert(Vector{Cint}, qr), convert(Vector{Cint}, qc), qv, rel, rhs)

function add_qconstr!(model::Model, lind::Vector{Cint}, lval::Vector{Float64}, qr::Vector{Cint}, qc::Vector{Cint}, qv::Vector{Float64}, rel::Char, rhs::Float64)
    qnnz = length(qr)
    if !(qnnz == length(qc) == length(qv))
        throw(ArgumentError("Inconsistent dimensions."))
    end

    lnnz = length(lind)
    if lnnz != length(lval)
        throw(ArgumentError("Inconsistent dimensions."))
    end
    
    if qnnz > 0
        ret = @grb_ccall(addqconstr, Cint, (
            Ptr{Void},    # model
            Cint,         # lnnz
            Ptr{Cint},    # lind
            Ptr{Float64}, # lval
            Cint,         # qnnz
            Ptr{Cint},    # qrow
            Ptr{Cint},    # qcol
            Ptr{Float64}, # qval
            Cchar,        # sense
            Float64,      # rhs
            Ptr{Uint8}    # name
            ), 
            model, lnnz, lind-1, lval, qnnz, qr-1, qc-1, qv, rel, rhs, C_NULL)
            
        if ret != 0
            throw(GurobiError(model.env, ret))
        end 
    end
    nothing
end


# add_rangeconstr

function add_rangeconstr!(model::Model, inds::Vector{Cint}, coeffs::Vector{Float64}, lower::Float64, upper::Float64)
   inds = inds - 1 # Zero-based indexing
   if !isempty(inds)
        ret = @grb_ccall(addrangeconstr, Cint, (
            Ptr{Void},    # model
            Cint,         # numnz
            Ptr{Cint},    # cind
            Ptr{Float64}, # cvals
            Float64,      # lower
			Float64,	  # upper
            Ptr{Uint8}    # name
            ),
            model, length(inds), inds, coeffs, lower, upper, C_NULL)
        if ret != 0
            throw(GurobiError(model.env, ret))
        end
    end
    nothing
end

function add_rangeconstrs!(model::Model, cbegins::Vector{Cint}, inds::Vector{Cint}, coeffs::Vector{Float64}, lower::Vector{Float64}, upper::Vector{Float64})
        
    m = length(cbegins)
    nnz = length(inds)
    
    if !(m == length(lower) == length(upper) && nnz == length(coeffs))
        throw(ArgumentError("Incompatible dimensions."))
    end 
        
    if m > 0 && nnz > 0
        ret = @grb_ccall(addrangeconstrs, Cint, (
            Ptr{Void},    # model
            Cint,         # num constraints
            Cint,         # num non-zeros
            Ptr{Cint},    # cbeg
            Ptr{Cint},    # cind
            Ptr{Float64}, # cval
            Ptr{Float64}, # lower
            Ptr{Float64}, # upper
            Ptr{Uint8}    # names
            ), 
            model, m, nnz, cbegins - 1, inds - 1, coeffs, 
            lower, upper, C_NULL)
        
        if ret != 0
            throw(GurobiError(model.env, ret))
        end
    end
    nothing
end


function _add_rangeconstrs_t!(model::Model, At::SparseMatrixCSC{Float64}, lower::Vector{Float64}, upper::Vector{Float64})
    cbeg = convert(Vector{Cint}, At.colptr[1:At.n])
    cind = convert(Vector{Cint}, At.rowval)
    add_rangeconstrs!(model, cbeg, cind, At.nzval, lower, upper)
end


function add_rangeconstrs!(model::Model, A::Matrix, lower::Vector{Float64}, upper::Vector{Float64})
    A = convert(Matrix{Float64},A)
    n::Int = num_vars(model)
    m::Int = size(A, 1)
    if !(m == length(upper) && n == size(A, 2))
        throw(ArgumentError("Incompatible dimensions."))
    end 
    
    At = sparse(transpose(A))  # each column of At now is a constraint
    _add_rangeconstrs_t!(model, At, lower, upper)
end

function add_rangeconstrs!(model::Model, A::SparseMatrixCSC, lower::Vector{Float64}, upper::Vector{Float64})
    A = convert(SparseMatrixCSC{Float64,Cint},A)
    n::Int = num_vars(model)
    m::Int = size(A, 1)
    if !(m == length(upper) && n == size(A, 2))
        throw(ArgumentError("Incompatible dimensions."))
    end 
    
    At = transpose(A)  # each column of At now is a constraint
    _add_rangeconstrs_t!(model, At, lower, upper)
end


#################################################
#
#  LP construction
#
#  minimize/maximize f'x
#
#   s.t.    A x <= b
#           Aeq x = beq
#           lb <= x <= ub
#
#################################################

typealias ConstrMat Union(Matrix{Float64}, SparseMatrixCSC{Float64})

function _add_constrs_t!(model::Model, At::SparseMatrixCSC{Float64}, rel::Char, b::Vector{Float64})    
    cbeg = convert(Vector{Cint}, At.colptr[1:At.n])
    cind = convert(Vector{Cint}, At.rowval)
    add_constrs!(model, cbeg, cind, At.nzval, rel, b)
end

function add_constrs!(model::Model, A::Matrix{Float64}, rel::Char, b::Vector{Float64})
    n::Int = num_vars(model)
    m::Int = size(A, 1)
    if !(m == length(b) && n == size(A, 2))
        throw(ArgumentError("Incompatible dimensions."))
    end 
    
    At = sparse(transpose(A))  # each column of At now is a constraint
    _add_constrs_t!(model, At, rel, b) 
end

function add_constrs!(model::Model, A::SparseMatrixCSC{Float64}, rel::Char, b::Vector{Float64})
    n::Int = num_vars(model)
    m::Int = size(A, 1)
    if !(m == length(b) && n == size(A, 2))
        throw(ArgumentError("Incompatible dimensions."))
    end 
    
    At = transpose(A)  # each column of At now is a constraint
    _add_constrs_t!(model, At, rel, b)
end


#################################################
#
#  QP construction
#
#  minimize (1/2) x'Hx + f'x
#
#   s.t.    A x <= b
#           Aeq x = beq
#           lb <= x <= ub
#
#################################################

function _add_qpterms!(model, H::Vector{Float64})  # H stores only the diagonal element
    n = num_vars(model)
    if n != length(H)
        throw(ArgumentError("Incompatible dimensions."))
    end
    qr = convert(Vector{Cint}, [0:n-1])
    qc = copy(qr)
    add_qpterms!(model, qr, qc, qv)
end

function _add_qpterms!(model, H::Float64)  # all diagonal elements are H
    n = num_vars(model)
    qr = convert(Vector{Cint}, [0:n-1])
    qc = copy(qr)
    qv = fill(H, n)
    add_qpterms!(model, qr, qc, qv)
end

function _add_qpterms!(model, H::SparseMatrixCSC{Float64}) # H must be symmetric
    n = num_vars(model)
    if !(H.m == n && H.n == n)
        throw(ArgumentError("H must be a symmetric matrix."))
    end
    
    nnz_h = nnz(H)
    qr = Array(Cint, nnz_h)
    qc = Array(Cint, nnz_h)
    qv = Array(Float64, nnz_h)
    k::Int = 0
    
    colptr::Vector{Int} = H.colptr
    rowval::Vector{Cint} = convert(Vector{Cint}, H.rowval)
    nzval::Vector{Float64} = convert(Vector{Float64}, H.nzval)
    
    for i = 1 : n
        qi::Cint = convert(Cint, i)
        for j = colptr[i]:(colptr[i+1]-1)
            qj::Cint = rowval[j] 
            
            if qi < qj
                k += 1
                qr[k] = qi
                qc[k] = qj
                qv[k] = nzval[j]
            elseif qi == qj
                k += 1
                qr[k] = qi
                qc[k] = qj
                qv[k] = nzval[j] * 0.5
            end
        end
    end
    
    add_qpterms!(model, qr[1:k], qc[1:k], qv[1:k])
end

function _add_qpterms!(model, H::Matrix{Float64}) # H must be symmetric
    n = num_vars(model)
    if !(size(H) == (n, n))
        throw(ArgumentError("H must be a symmetric matrix."))
    end
    
    nmax = int(n * (n + 1) / 2)
    qr = Array(Cint, nmax)
    qc = Array(Cint, nmax)
    qv = Array(Float64, nmax)
    k::Int = 0
    
    for i = 1 : n
        qi::Cint = convert(Cint, i)
        
        v::Float64 = H[i,i]
        if v != 0.
            k += 1
            qr[k] = qi
            qc[k] = qi
            qv[k] = v * 0.5
        end
        
        for j = i+1 : n
            v = H[j, i]
            if v != 0.
                k += 1
                qr[k] = qi
                qc[k] = convert(Cint, j)
                qv[k] = v
            end
        end
    end
        
    add_qpterms!(model, qr[1:k], qc[1:k], qv[1:k])
end


function qp_model(env::Env, name::ASCIIString, 
    H::Union(Vector{Float64}, Matrix{Float64}, SparseMatrixCSC{Float64}, Float64), 
    f::Vector{Float64}, 
    A::Union(ConstrMat, Nothing), 
    b::Union(Vector{Float64}, Nothing), 
    Aeq::Union(ConstrMat, Nothing), 
    beq::Union(Vector{Float64}, Nothing), 
    lb::Bounds, ub::Bounds)
    
    # create model
    model = gurobi_model(env, name)
    
    # add variables
    add_cvars!(model, f, lb, ub)
    update_model!(model)
    
    # add qpterms
    
    _add_qpterms!(model, H)
    
    # add constraints
    if A != nothing && b != nothing
        add_constrs!(model, A, '<', b)
    end
    
    if Aeq != nothing && beq != nothing
        add_constrs!(model, Aeq, '=', beq)
    end
    update_model!(model)
    
    model
end

qp_model(env::Env, name, H, f, A, b, Aeq, beq) = qp_model(env, name, H, f, A, b, Aeq, beq, nothing, nothing)
qp_model(env::Env, name, H, f, A, b) = qp_model(env, name, H, f, A, b, nothing, nothing, nothing, nothing)

