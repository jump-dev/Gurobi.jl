# Gurobi model attributes

############################################
#
#   Low level attribute getter/setters
#
############################################

function get_intattr(model::Model, name::ASCIIString)
    a = Array(Cint, 1)
    ret = @grb_ccall(getintattr, Cint, 
        (Ptr{Void}, Ptr{Uint8}, Ptr{Cint}),
        model, name, a);
    if ret != 0
        throw(GurobiError(model.env, ret))
    end
    convert(Int, a[1])
end

function get_dblattr(model::Model, name::ASCIIString)
    a = Array(Float64, 1)
    ret = @grb_ccall(getdblattr, Cint, 
        (Ptr{Void}, Ptr{Uint8}, Ptr{Float64}),
        model, name, a);
    if ret != 0
        throw(GurobiError(model.env, ret))
    end
    a[1]::Float64
end

function get_strattr(model::Model, name::ASCIIString)
    a = Array(Ptr{Uint8}, 1)
    ret = @grb_ccall(getstrattr, Cint, 
        (Ptr{Void}, Ptr{Uint8}, Ptr{Ptr{Uint8}}), 
        model, name, a)
    if ret != 0
        throw(GurobiError(model.env, ret))
    end
    bytestring(a[1])
end

# Note: in attrarray API, the start argument is one-based (following Julia convention)

function get_intattrarray!(r::Array{Cint}, model::Model, name::ASCIIString, start::Integer)
    ret = @grb_ccall(getintattrarray, Cint, 
        (Ptr{Void}, Ptr{Uint8}, Cint, Cint, Ptr{Cint}), 
        model, name, start - 1, length(r), r)
    if ret != 0
        throw(GurobiError(model.env, ret))
    end
    r    
end

function get_intattrarray(model::Model, name::ASCIIString, start::Integer, len::Integer)
    get_intattrarray!(Array(Cint, len), model, name, start)
end

function get_dblattrarray!(r::Array{Float64}, model::Model, name::ASCIIString, start::Integer)
    ret = @grb_ccall(getdblattrarray, Cint, 
        (Ptr{Void}, Ptr{Uint8}, Cint, Cint, Ptr{Float64}), 
        model, name, start - 1, length(r), r)
    if ret != 0
        throw(GurobiError(model.env, ret))
    end
    r    
end

function get_dblattrarray(model::Model, name::ASCIIString, start::Integer, len::Integer)
    get_dblattrarray!(Array(Float64, len), model, name, start)
end

function get_charattrarray!(r::Array{Cchar}, model::Model, name::ASCIIString, start::Integer)
    ret = @grb_ccall(getcharattrarray, Cint, 
        (Ptr{Void}, Ptr{Uint8}, Cint, Cint, Ptr{Cchar}), 
        model, name, start - 1, length(r), r)
    if ret != 0
        throw(GurobiError(model.env, ret))
    end
    char(r)
end

function get_charattrarray(model::Model, name::ASCIIString, start::Integer, len::Integer)
    get_charattrarray!(Array(Cchar, len), model, name, start)
end

# setters

function set_intattr!(model::Model, name::ASCIIString, v::Integer)
    ret = @grb_ccall(setintattr, Cint, 
        (Ptr{Void}, Ptr{Uint8}, Cint), model, name, v)
    if ret != 0
        throw(GurobiError(model.env, ret))
    end
    nothing
end

function set_dblattr!(model::Model, name::ASCIIString, v::Real)
    ret = @grb_ccall(setdblattr, Cint, 
        (Ptr{Void}, Ptr{Uint8}, Float64), model, name, v)
    if ret != 0
        throw(GurobiError(model.env, ret))
    end
    nothing
end

function set_strattr!(model::Model, name::ASCIIString, v::ASCIIString)
    ret = @grb_ccall(setstrattr, Cint, 
        (Ptr{Void}, Ptr{Uint8}, Ptr{Uint8}), model, name, v)
    if ret != 0
        throw(GurobiError(model.env, ret))
    end
    nothing
end

# reminder: start is one-based

function set_intattrarray!(model::Model, name::ASCIIString, start::Integer, len::Integer, values::Vector)
    values = convert(Vector{Cint},values)
    ret = @grb_ccall(setintattrarray, Cint, 
        (Ptr{Void}, Ptr{Uint8}, Cint, Cint, Ptr{Cint}), model, name, start-1, len, ivec(values))
    if ret != 0
        throw(GurobiError(model.env, ret))
    end
    nothing
end

function set_dblattrarray!(model::Model, name::ASCIIString, start::Integer, len::Integer, values::Vector)
    values = convert(Vector{Float64},values)
    ret = @grb_ccall(setdblattrarray, Cint, 
        (Ptr{Void}, Ptr{Uint8}, Cint, Cint, Ptr{Float64}), model, name, start-1, len, fvec(values))
    if ret != 0
        throw(GurobiError(model.env, ret))
    end
    nothing
end

function set_charattrarray!(model::Model, name::ASCIIString, start::Integer, len::Integer, values::Vector)
    values = convert(Vector{Cchar},values)
    ret = @grb_ccall(setcharattrarray, Cint, 
        (Ptr{Void}, Ptr{Uint8}, Cint, Cint, Ptr{Cchar}), model, name, start-1, len, cvec(values))
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

