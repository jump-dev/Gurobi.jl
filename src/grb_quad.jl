# Quadratic terms & constraints
#

function add_qpterms!(model::Model, qr::IVec, qc::IVec, qv::FVec)
    nnz = length(qr)
    (nnz == length(qc) == length(qv)) || error("Inconsistent argument dimensions.")

    if nnz > 0
        ret = @grb_ccall(addqpterms, Cint, (
            Ptr{Cvoid},    # model
            Cint,         # nnz
            Ptr{Cint},    # qrow
            Ptr{Cint},    # qcol
            Ptr{Float64}, # qval
            ),
            model, nnz, qr.-Cint(1), qc.-Cint(1), qv)

        if ret != 0
            throw(GurobiError(model.env, ret))
        end
    end
    nothing
end

function add_qpterms!(model::Model, qr::Vector, qc::Vector, qv::Vector)
    add_qpterms!(model, ivec(qr), ivec(qc), fvec(qv))
end


function add_qpterms!(model, H::SparseMatrixCSC{Float64}) # H must be symmetric
    n = num_vars(model)
    (H.m == n && H.n == n) || error("H must be an n-by-n symmetric matrix.")

    nnz_h = nnz(H)
    qr = Array{Cint}(undef, nnz_h)
    qc = Array{Cint}(undef, nnz_h)
    qv = Array{Float64}(undef, nnz_h)
    k = 0

    colptr::Vector{Int} = H.colptr
    nzval::Vector{Float64} = H.nzval

    for i = 1 : n
        qi::Cint = convert(Cint, i)
        for j = colptr[i]:(colptr[i+1]-1)
            qj = convert(Cint, H.rowval[j])

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
    size(H) == (n, n) || error("H must be an n-by-n symmetric matrix.")

    nmax = round(Int,n * (n + 1) / 2)
    qr = Array{Cint}(undef, nmax)
    qc = Array{Cint}(undef, nmax)
    qv = Array{Float64}(undef, nmax)
    k::Int = 0

    for i = 1 : n
        qi = convert(Cint, i)

        v = H[i,i]
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

function add_diag_qpterms!(model, H::Vector)  # H stores only the diagonal element
    n = num_vars(model)
    n == length(H) || error("Incompatible dimensions.")
    q = [convert(Cint,1):convert(Cint,n)]
    add_qpterms!(model, q, q, fvec(h))
end

function add_diag_qpterms!(model, hv::Real)  # all diagonal elements are H
    n = num_vars(model)
    q = [convert(Cint,1):convert(Cint,n)]
    add_qpterms!(model, q, q, fill(float64(hv), n))
end

function delq!(model::Model)
    ret = @grb_ccall(delq, Cint, (
        Ptr{Cvoid},    # model
        ),
        model)

    if ret != 0
        throw(GurobiError(model.env, ret))
    end
end

function getq(model::Model)
    nz = get_intattr(model, "NumQNZs")
    rowidx = Array{Cint}(undef, nz)
    colidx = Array{Cint}(undef, nz)
    val = Array{Float64}(undef, nz)
    nzout = Ref{Cint}()

    ret = @grb_ccall(getq, Cint, (
        Ptr{Cvoid},  # model
        Ptr{Cint},  # numqnzP
        Ptr{Cint},  # qrow
        Ptr{Cint},  # qcol
        Ptr{Float64}# qval
        ),
        model,nzout,rowidx,colidx,val)

    if ret != 0
        throw(GurobiError(model.env, ret))
    end

    return rowidx, colidx, val
end

# add_qconstr!

function add_qconstr!(model::Model, lind::IVec, lval::FVec, qr::IVec, qc::IVec, qv::FVec, rel::Cchar, rhs::Float64)
    qnnz = length(qr)
    qnnz == length(qc) == length(qv) || error("Inconsistent argument dimensions.")

    lnnz = length(lind)
    lnnz == length(lval) || error("Inconsistent argument dimensions.")

    if qnnz > 0 || lnnz > 0
        ret = @grb_ccall(addqconstr, Cint, (
            Ptr{Cvoid},    # model
            Cint,         # lnnz
            Ptr{Cint},    # lind
            Ptr{Float64}, # lval
            Cint,         # qnnz
            Ptr{Cint},    # qrow
            Ptr{Cint},    # qcol
            Ptr{Float64}, # qval
            Cchar,        # sense
            Float64,      # rhs
            Ptr{UInt8}    # name
            ),
            model, lnnz, lind.-Cint(1), lval, qnnz, qr.-Cint(1), qc.-Cint(1), qv, rel, rhs, C_NULL)

        if ret != 0
            throw(GurobiError(model.env, ret))
        end
    end
    nothing
end

function add_qconstr!(model::Model, lind::Vector, lval::Vector, qr::Vector, qc::Vector,
    qv::Vector{Float64}, rel::GChars, rhs::Real)

    add_qconstr!(model, ivec(lind), fvec(lval), ivec(qr), ivec(qc), fvec(qv), cchar(rel), Float64(rhs))
end

function delqconstrs!(model::Model, indices::Vector{Int})
    ret = @grb_ccall(delqconstrs, Cint, (Ptr{Cvoid}, Cint, Ptr{Cint}), model, length(indices), Cint.(indices.-1))
    if ret != 0
        throw(GurobiError(model.env, ret))
    end
end

function getqconstr(model::Model, index::Int)
    # From the Gurobi documentation:
    # Typical usage is to call this routine twice. In the first call, you
    # specify the requested quadratic constraint, with NULL values for the array
    # arguments.
    affine_nnz    = Cint[0]
    quadratic_nnz = Cint[0]
    ret = @grb_ccall(getqconstr, Cint, (Ptr{Cvoid}, Cint,  # model, constraint
        Ptr{Cint}, Ptr{Cint}, Ptr{Float64},                  # affine terms
        Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Float64}        # quadratic terms
        ),
        model, Cint(index-1),
        affine_nnz, C_NULL, C_NULL,
        quadratic_nnz, C_NULL, C_NULL, C_NULL
    )
    if ret != 0
        throw(GurobiError(model.env, ret))
    end
    affine_cols     = fill(Cint(0), affine_nnz[1])
    affine_coefs    = fill(0.0, affine_nnz[1])
    quadratic_rows  = fill(Cint(0), quadratic_nnz[1])
    quadratic_cols  = fill(Cint(0), quadratic_nnz[1])
    quadratic_coefs = fill(0.0, quadratic_nnz[1])

    ret = @grb_ccall(getqconstr, Cint, (Ptr{Cvoid}, Cint,  # model, constraint
        Ptr{Cint}, Ptr{Cint}, Ptr{Float64},                  # affine terms
        Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Float64}        # quadratic terms
        ),
        model, Cint(index-1),
        affine_nnz, affine_cols, affine_coefs,
        quadratic_nnz, quadratic_rows, quadratic_cols, quadratic_coefs
    )
    if ret != 0
        throw(GurobiError(model.env, ret))
    end
    return affine_cols, affine_coefs, quadratic_rows, quadratic_cols, quadratic_coefs
end
