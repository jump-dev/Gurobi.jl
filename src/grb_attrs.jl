# Gurobi model attributes

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

