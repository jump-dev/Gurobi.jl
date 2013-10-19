# Add constraints

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


# add_rangeconstr! & add_rangeconstrs!

function add_rangeconstr!(model::Model, inds::Vector{Cint}, coeffs::Vector{Float64}, lower::Float64, upper::Float64)
   inds = inds - 1 # Zero-based indexing
   if !isempty(inds)
        ret = @grb_ccall(addrangeconstr, Cint, (
            Ptr{Void},    # model
            Cint,         # numnz
            Ptr{Cint},    # cind
            Ptr{Float64}, # cvals
            Float64,      # lower
            Float64,      # upper
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

