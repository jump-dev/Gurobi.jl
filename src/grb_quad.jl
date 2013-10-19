# Quadratic terms & constraints
#

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

function add_qpterms!(model, H::SparseMatrixCSC{Float64}) # H must be symmetric
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

function add_qpterms!(model, H::Matrix{Float64}) # H must be symmetric
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

function add_diag_qpterms!(model, H::Vector{Float64})  # H stores only the diagonal element
    n = num_vars(model)
    if n != length(H)
        throw(ArgumentError("Incompatible dimensions."))
    end
    qr = convert(Vector{Cint}, [0:n-1])
    qc = copy(qr)
    add_qpterms!(model, qr, qc, qv)
end

function add_diag_qpterms!(model, H::Float64)  # all diagonal elements are H
    n = num_vars(model)
    qr = convert(Vector{Cint}, [0:n-1])
    qc = copy(qr)
    qv = fill(H, n)
    add_qpterms!(model, qr, qc, qv)
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