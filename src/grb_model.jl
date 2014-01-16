# Gurobi model

#################################################
#
#  model type & constructors
#
#################################################

type Model
    env::Env
    ptr_model::Ptr{Void}
    callback::Any
    
    function Model(env::Env, p::Ptr{Void})
        model = new(env, p, nothing)
        finalizer(model, free_model)
        model
    end
end

function Model(env::Env, name::ASCIIString)

    @assert is_valid(env)
    
    a = Array(Ptr{Void}, 1)
    ret = @grb_ccall(newmodel, Cint, (
        Ptr{Void},  # env
        Ptr{Ptr{Void}},  # pointer to model pointer
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


function Model(env::Env, name::ASCIIString, sense::Symbol)
    model = Model(env, name)
    if sense != :minimize
        set_sense!(model, sense)
    end
    model 
end


#################################################
#
#  model manipulation
#
#################################################

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

function read_model(model::Model, filename::ASCIIString)
    @assert is_valid(model.env)
    a = Array(Ptr{Void}, 1)
    ret = @grb_ccall(readmodel, Cint, 
        (Ptr{Void}, Ptr{Uint8}, Ptr{Ptr{Void}}), 
        model.env, filename, a)
    if ret != 0
        throw(GurobiError(model.env, ret))
    end
    model.ptr_model = a[1]
    nothing
end

function write_model(model::Model, filename::ASCIIString)
    ret = @grb_ccall(write, Cint, (Ptr{Void}, Ptr{Uint8}), 
        model.ptr_model, filename)
    if ret != 0
        throw(GurobiError(model.env, ret))
    end
    nothing
end

