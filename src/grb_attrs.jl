# Gurobi model attributes

############################################
#
#   Low level attribute getter/setters
#
############################################

function get_intattr(model::Model, name::String)
    @assert isascii(name)
    a = Ref{Cint}()
    ret = @grb_ccall(getintattr, Cint,
        (Ptr{Cvoid}, Ptr{UInt8}, Ref{Cint}),
        model, name, a);
    if ret != 0
        throw(GurobiError(model.env, ret))
    end
    Int(a[])
end

function get_dblattr(model::Model, name::String)
    @assert isascii(name)
    a = Ref{Float64}()
    ret = @grb_ccall(getdblattr, Cint,
        (Ptr{Cvoid}, Ptr{UInt8}, Ref{Float64}),
        model, name, a);
    if ret != 0
        throw(GurobiError(model.env, ret))
    end
    a[]
end

function get_strattr(model::Model, name::String)
    @assert isascii(name)
    a = Ref{Ptr{UInt8}}()
    ret = @grb_ccall(getstrattr, Cint,
        (Ptr{Cvoid}, Ptr{UInt8}, Ref{Ptr{UInt8}}),
        model, name, a)
    if ret != 0
        throw(GurobiError(model.env, ret))
    end
    unsafe_string(a[])
end

# array element

function get_intattrelement(model::Gurobi.Model, name::String, element::Int)
    @assert isascii(name)
    a = Ref{Cint}()
    ret = @grb_ccall(getintattrelement, Cint,
        (Ptr{Cvoid}, Ptr{UInt8}, Cint, Ref{Cint}),
        model, name, element - 1, a
    )
    if ret != 0
        throw(GurobiError(model.env, ret))
    end
    Int(a[])
end

function get_dblattrelement(model::Gurobi.Model, name::String, element::Int)
    @assert isascii(name)
    a = Ref{Float64}()
    ret = @grb_ccall(getdblattrelement, Cint,
        (Ptr{Cvoid}, Ptr{UInt8}, Cint, Ref{Float64}),
        model, name, element - 1, a
    )
    if ret != 0
        throw(GurobiError(model.env, ret))
    end
    a[]
end

function get_charattrelement(model::Gurobi.Model, name::String, element::Int)
    @assert isascii(name)
    a = Ref{Cchar}()
    ret = @grb_ccall(getcharattrelement, Cint,
        (Ptr{Cvoid}, Ptr{UInt8}, Cint, Ref{Cchar}),
        model, name, element - 1, a
    )
    if ret != 0
        throw(GurobiError(model.env, ret))
    end
    Char(a[])
end

function get_strattrelement(model::Gurobi.Model, name::String, element::Int)
    @assert isascii(name)
    a = Ref{Ptr{UInt8}}()
    ret = @grb_ccall(getstrattrelement, Cint,
        (Ptr{Cvoid}, Ptr{UInt8}, Cint, Ref{Ptr{UInt8}}),
        model, name, element - 1, a
    )
    if ret != 0
        throw(GurobiError(model.env, ret))
    end
    return unsafe_string(a[])
end

# Note: in attrarray API, the start argument is one-based (following Julia convention)

function get_intattrarray!(r::Array{Cint}, model::Model, name::String, start::Integer)
    @assert isascii(name)
    ret = @grb_ccall(getintattrarray, Cint,
        (Ptr{Cvoid}, Ptr{UInt8}, Cint, Cint, Ptr{Cint}),
        model, name, start - 1, length(r), r)
    if ret != 0
        throw(GurobiError(model.env, ret))
    end
    r
end

function get_intattrlist!(r::Array{Cint}, model::Model, name::String, inds::Vector{I}) where I<:Integer
    @assert isascii(name)
    @assert length(r) == length(inds)
    ret = @grb_ccall(getintattrlist, Cint,
        (Ptr{Cvoid}, Ptr{UInt8}, Cint, Ptr{Cint}, Ptr{Cint}),
        model, name, Cint(length(inds)), ivec(inds .- 1), r)
    if ret != 0
        throw(GurobiError(model.env, ret))
    end
    nothing
end

function get_intattrlist(model::Model, name::String, inds::Vector{I}) where I<:Integer
    r = Array{Cint}(undef, length(inds))
    get_intattrlist!(r::Array{Cint}, model::Model, name::String, inds)
    return r
end

function get_intattrarray(model::Model, name::String, start::Integer, len::Integer)
    @assert isascii(name)
    get_intattrarray!(Array{Cint}(len), model, name, start)
end

function get_dblattrarray!(r::Array{Float64}, model::Model, name::String, start::Integer)
    @assert isascii(name)
    ret = @grb_ccall(getdblattrarray, Cint,
        (Ptr{Cvoid}, Ptr{UInt8}, Cint, Cint, Ptr{Float64}),
        model, name, start - 1, length(r), r)
    if ret != 0
        throw(GurobiError(model.env, ret))
    end
    r
end

function get_dblattrlist!(r::Array{Float64}, model::Model, name::String, inds::Vector{I}) where I<:Integer
    @assert isascii(name)
    @assert length(r) == length(inds)
    ret = @grb_ccall(getdblattrlist, Cint,
        (Ptr{Cvoid}, Ptr{UInt8}, Cint, Ptr{Cint}, Ptr{Float64}),
        model, name, Cint(length(inds)), ivec(inds .- 1), r)
    if ret != 0
        throw(GurobiError(model.env, ret))
    end
    nothing
end

function get_dblattrlist(model::Model, name::String, inds::Vector{I}) where I<:Integer
    r = Array{Float64}(undef, length(inds))
    get_dblattrlist!(r::Array{Float64}, model::Model, name::String, inds)
    return r
end

function get_dblattrarray(model::Model, name::String, start::Integer, len::Integer)
    @assert isascii(name)
    get_dblattrarray!(Array{Float64}(undef, len), model, name, start)
end

function get_charattrarray!(r::Array{Cchar}, model::Model, name::String, start::Integer)
    @assert isascii(name)
    ret = @grb_ccall(getcharattrarray, Cint,
        (Ptr{Cvoid}, Ptr{UInt8}, Cint, Cint, Ptr{Cchar}),
        model, name, start - 1, length(r), r)
    if ret != 0
        throw(GurobiError(model.env, ret))
    end
    map(Char,r)
end

function get_charattrlist!(r::Array{Cchar}, model::Model, name::String, inds::Vector{I}) where I<:Integer
    @assert isascii(name)
    @assert length(r) == length(inds)
    ret = @grb_ccall(getcharattrlist, Cint,
        (Ptr{Cvoid}, Ptr{UInt8}, Cint, Ptr{Cint}, Ptr{Cchar}),
        model, name, Cint(length(inds)), ivec(inds .- 1), r)
    if ret != 0
        throw(GurobiError(model.env, ret))
    end
    nothing
end

function get_charattrlist(model::Model, name::String, inds::Vector{I}) where I<:Integer
    r = Array{Cchar}(undef, length(inds))
    get_charattrlist!(r::Array{Cchar}, model::Model, name::String, inds)
    return r
end

function get_charattrarray(model::Model, name::String, start::Integer, len::Integer)
    @assert isascii(name)
    get_charattrarray!(Array{Cchar}(undef, len), model, name, start)
end

function get_strattrarray(model::Model, name::String, start::Integer, len::Integer)
    @assert isascii(name)
    get_strattrarray!(Array{Ptr{UInt8}}(undef, len), model, name, start)
end

function get_strattrarray!(r::Array{Ptr{UInt8}}, model::Model, name::String, start::Integer)
    @assert isascii(name)
    ret = @grb_ccall(getstrattrarray, Cint,
        (Ptr{Cvoid}, Ptr{UInt8}, Cint, Cint, Ptr{Ptr{UInt8}}),
        model, name, start - 1, length(r), r)
    if ret != 0
        throw(GurobiError(model.env, ret))
    end
    map(unsafe_string, r)
end

# setters

function set_intattr!(model::Model, name::String, v::Integer)
    @assert isascii(name)
    ret = @grb_ccall(setintattr, Cint,
        (Ptr{Cvoid}, Ptr{UInt8}, Cint), model, name, v)
    if ret != 0
        throw(GurobiError(model.env, ret))
    end
    nothing
end

function set_dblattr!(model::Model, name::String, v::Real)
    @assert isascii(name)
    ret = @grb_ccall(setdblattr, Cint,
        (Ptr{Cvoid}, Ptr{UInt8}, Float64), model, name, v)
    if ret != 0
        throw(GurobiError(model.env, ret))
    end
    nothing
end

function set_strattr!(model::Model, name::String, v::String)
    @assert isascii(name)
    ret = @grb_ccall(setstrattr, Cint,
        (Ptr{Cvoid}, Ptr{UInt8}, Ptr{UInt8}), model, name, v)
    if ret != 0
        throw(GurobiError(model.env, ret))
    end
    nothing
end

function set_strattrelement!(model::Model, name::String, el::Int, v::String)
    @assert isascii(name)
    ret = @grb_ccall(setstrattrelement, Cint,
        (Ptr{Cvoid}, Ptr{UInt8}, Cint, Ptr{UInt8}), model, name, el - 1, v)
    if ret != 0
        throw(GurobiError(model.env, ret))
    end
    return
end

# array element

function set_intattrelement!(model::Gurobi.Model, name::String, element::Int, v::Integer)
    @assert isascii(name)
    ret = @grb_ccall(setintattrelement, Cint,
        (Ptr{Cvoid}, Ptr{UInt8}, Cint, Cint),
        model, name, element - 1, v
    )
    if ret != 0
        throw(GurobiError(model.env, ret))
    end
    nothing
end

function set_dblattrelement!(model::Gurobi.Model, name::String, element::Int, v::Real)
    @assert isascii(name)
    ret = @grb_ccall(setdblattrelement, Cint,
        (Ptr{Cvoid}, Ptr{UInt8}, Cint, Float64),
        model, name, element - 1, v
    )
    if ret != 0
        throw(GurobiError(model.env, ret))
    end
    nothing
end

function set_charattrelement!(model::Gurobi.Model, name::String, element::Int, v::Char)
    @assert isascii(name)
    ret = @grb_ccall(setcharattrelement, Cint,
        (Ptr{Cvoid}, Ptr{UInt8}, Cint, Cchar),
        model, name, element - 1, v
    )
    if ret != 0
        throw(GurobiError(model.env, ret))
    end
    nothing
end

# reminder: start is one-based

function set_intattrarray!(model::Model, name::String, start::Integer, len::Integer, values::Vector)
    @assert isascii(name)
    values = convert(Vector{Cint},values)
    ret = @grb_ccall(setintattrarray, Cint,
        (Ptr{Cvoid}, Ptr{UInt8}, Cint, Cint, Ptr{Cint}), model, name, start-1, len, ivec(values))
    if ret != 0
        throw(GurobiError(model.env, ret))
    end
    nothing
end

function set_intattrlist!(model::Model, name::String, inds::Vector, values::Vector)
    @assert isascii(name)
    @assert length(inds) == length(values)
    ret = @grb_ccall(setintattrlist, Cint,
        (Ptr{Cvoid}, Ptr{UInt8}, Cint, Ptr{Cint}, Ptr{Cint}), model, name, Cint(length(inds)), ivec(inds.-1), ivec(values))
    if ret != 0
        throw(GurobiError(model.env, ret))
    end
    nothing
end

function set_dblattrarray!(model::Model, name::String, start::Integer, len::Integer, values::Vector)
    @assert isascii(name)
    values = convert(Vector{Float64},values)
    ret = @grb_ccall(setdblattrarray, Cint,
        (Ptr{Cvoid}, Ptr{UInt8}, Cint, Cint, Ptr{Float64}), model, name, start-1, len, fvec(values))
    if ret != 0
        throw(GurobiError(model.env, ret))
    end
    nothing
end

function set_dblattrlist!(model::Model, name::String, inds::Vector, values::Vector)
    @assert isascii(name)
    @assert length(inds) == length(values)
    ret = @grb_ccall(setdblattrlist, Cint,
        (Ptr{Cvoid}, Ptr{UInt8}, Cint, Ptr{Cint}, Ptr{Float64}), model, name, Cint(length(inds)), ivec(inds.-1), fvec(values))
    if ret != 0
        throw(GurobiError(model.env, ret))
    end
    nothing
end

function set_charattrarray!(model::Model, name::String, start::Integer, len::Integer, values::Vector)
    @assert isascii(name)
    values = convert(Vector{Cchar},values)
    ret = @grb_ccall(setcharattrarray, Cint,
        (Ptr{Cvoid}, Ptr{UInt8}, Cint, Cint, Ptr{Cchar}), model, name, start-1, len, cvec(values))
    if ret != 0
        throw(GurobiError(model.env, ret))
    end
    nothing
end

function set_charattrlist!(model::Model, name::String, inds::Vector, values::Vector)
    @assert isascii(name)
    @assert length(inds) == length(values)
    ret = @grb_ccall(setcharattrlist, Cint,
        (Ptr{Cvoid}, Ptr{UInt8}, Cint, Ptr{Cint}, Ptr{Cchar}), model, name, Cint(length(inds)), ivec(inds.-1), cvec(values))
    if ret != 0
        throw(GurobiError(model.env, ret))
    end
    nothing
end

############################################
#
#   Macros for array definition
#
############################################

macro grb_int_attr(fun, attrname)
    @eval $(fun)(model::Model) = get_intattr(model, $(attrname))
end

macro grb_dbl_attr(fun, attrname)
    @eval $(fun)(model::Model) = get_dblattr(model, $(attrname))
end

macro grb_str_attr(fun, attrname)
    @eval $(fun)(model::Model) = get_strattr(model, $(attrname))
end


############################################
#
#   Model attributes
#
############################################

# basic attributes

@grb_str_attr model_name  "ModelName"

@grb_int_attr num_vars     "NumVars"
@grb_int_attr num_constrs  "NumConstrs"
@grb_int_attr num_sos      "NumSOS"
@grb_int_attr num_qconstrs "NumQConstrs"
@grb_int_attr num_cnzs     "NumNZs"
@grb_int_attr num_qnzs     "NumQNZs"
@grb_int_attr num_qcnzs    "NumQCNZs"

@grb_int_attr num_intvars  "NumIntVars"
@grb_int_attr num_binvars  "NumBinVars"

# derived attribute functions

model_sense(model::Model) = get_intattr(model, "ModelSense") > 0 ? (:minimize) : (:maximize)

is_mip(model::Model) = get_intattr(model, "IsMIP") != 0
is_qp(model::Model)  = get_intattr(model, "IsQP") != 0
is_qcp(model::Model) = get_intattr(model, "IsQCP") != 0

function model_type(model::Model)
    is_qp(model)  ? (:QP)  :
    is_qcp(model) ? (:QCP) : (:LP)
end

function set_sense!(model::Model, sense::Symbol)
    v = sense == :maximize ? -1 :
        sense == :minimize ? 1 :
        throw(ArgumentError("Invalid model sense."))

    set_intattr!(model, "ModelSense", v)
end

# variable related attributes

lowerbounds(model::Model) = get_dblattrarray(model, "LB", 1, num_vars(model))
upperbounds(model::Model) = get_dblattrarray(model, "UB", 1, num_vars(model))
objcoeffs(model::Model) = get_dblattrarray(model, "Obj", 1, num_vars(model))

# note: this takes effect only after update_model! is called:
function set_objcoeffs!(model::Model, c::Vector)
    n = num_vars(model)
    length(c) == n || error("Inconsistent argument dimensions.")
    set_dblattrarray!(model, "Obj", 1, n, c)
end


############################################
#
#   The show method for model
#
#   - Based on attributes
#
############################################

function show(io::IO, model::Model)
    if model.ptr_model != C_NULL
        println(io, "Gurobi Model: $(model_name(model))")
        if is_mip(model)
            println(io, "    type   : $(model_type(model)) (MIP)")
        else
            println(io, "    type   : $(model_type(model))")
        end
        println(io, "    sense  : $(model_sense(model))")
        println(io, "    number of variables             = $(num_vars(model))")
        println(io, "    number of linear constraints    = $(num_constrs(model))")
        println(io, "    number of quadratic constraints = $(num_qconstrs(model))")
        println(io, "    number of sos constraints       = $(num_sos(model))")
        println(io, "    number of non-zero coeffs       = $(num_cnzs(model))")
        println(io, "    number of non-zero qp objective terms  = $(num_qnzs(model))")
        println(io, "    number of non-zero qp constraint terms = $(num_qcnzs(model))")
    else
        println(io, "Gurobi Model: NULL")
    end
end
