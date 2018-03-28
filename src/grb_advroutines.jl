type GRBsvec
    len::Cint
    ind::Ptr{Cint}
    val::Ptr{Cdouble}
end

function get_tableaurow(model::Model, cidx::Integer)
    len = num_vars(model)+num_constrs(model)
    grb_v = GRBsvec(0, pointer(Array{Cint}(len)) , pointer(Array{Cdouble}(len) ))
    ret = @grb_ccall(BinvRowi, Cint, (
        Ptr{Void},
        Cint,
        Ptr{GRBsvec}),
        model, convert(Cint, cidx-1), &grb_v)
    if ret != 0
        throw(GurobiError(model.env, ret))
    end

    idx = unsafe_wrap(Array, grb_v.ind, grb_v.len)
    val = unsafe_wrap(Array, grb_v.val, grb_v.len)
    spv = spzeros(len)
    spv[idx+1] = val
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

# TODO GRBBinvi
