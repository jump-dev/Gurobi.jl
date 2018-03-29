mutable struct GRBsvec
    len::Cint
    ind::Ptr{Cint}
    val::Ptr{Cdouble}
end

function get_tableaurow_grb(model::Model, cidx::Integer)
	len = num_vars(model)+num_constrs(model)
	idx = Array{Cint}(len)
	val = Array{Cdouble}(len)

    grb_v = GRBsvec(0, pointer(idx) , pointer(val))
    ret = @grb_ccall(BinvRowi, Cint, (
        Ptr{Void},
        Cint,
        Ref{GRBsvec}),
        model, convert(Cint, cidx-1), grb_v)
    if ret != 0
        throw(GurobiError(model.env, ret))
    end
	return idx[1:grb_v.len]+1, val[1:grb_v.len]
end

function get_tableaurow!(row::Vector{Float64}, model::Model, cidx::Integer)
	idx, val = get_tableaurow_grb(model, cidx)
    row[idx] = val
    return nothing
end

function get_tableaurow!(BinvA::Matrix{Float64}, model::Model, cidx::Integer)
    idx, val = get_tableaurow_grb(model, cidx)
    BinvA[cidx, idx] = val
    return nothing
end

function get_tableaurow(model::Model, cidx::Integer)
	len = num_vars(model)+num_constrs(model)
	spv = zeros(len)
	get_tableaurow!(spv, model, cidx)
	return spv
end

function get_basisidx(model::Model)
    v = Array{Cint}(num_constrs(model))
    ret = @grb_ccall(getBasisHead, Cint, (
        Ptr{Void},
        Ptr{Int}),
        model, v)
    if ret != 0
        throw(GurobiError(model.env, ret))
    end
    return v+1
end
