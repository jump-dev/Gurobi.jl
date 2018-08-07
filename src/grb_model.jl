# Gurobi model

#################################################
#
#  model type & constructors
#
#################################################

mutable struct Model
    env::Env
    ptr_model::Ptr{Cvoid}
    callback::Any
    finalize_env::Bool

    function Model(env::Env, p::Ptr{Cvoid}; finalize_env::Bool=false)
        model = new(env, p, nothing, finalize_env)
        if VERSION >= v"0.7-"
            finalizer(model) do m
                free_model(m)
                if m.finalize_env
                    free_env(m.env)
                end
            end
        else
            finalizer(model, m -> (free_model(m); if m.finalize_env; free_env(m.env); end))
        end
        model
    end
end

# If environment is tied to this model, the, we can safely finalize
# both at the same time, working around the Julia GC.
function Model(env::Env, name::String; finalize_env::Bool=false)

    @assert is_valid(env)
    @assert isascii(name)

    a = Ref{Ptr{Cvoid}}()
    ret = @grb_ccall(newmodel, Cint, (
        Ptr{Cvoid},  # env
        Ptr{Ptr{Cvoid}},  # pointer to model pointer
        Ptr{UInt8},  # name
        Cint,  # numvars
        Ptr{Float64}, # obj coeffs
        Ptr{Float64}, # lbounds
        Ptr{Float64}, # ubounds
        Ptr{UInt8},   # var types,
        Ptr{Ptr{UInt8}} # varnames
        ),
        env, a, name, 0,
        C_NULL, C_NULL, C_NULL, C_NULL, C_NULL
    )

    if ret != 0
        throw(GurobiError(env, ret))
    end

    Model(env, a[]; finalize_env=finalize_env)
end


function Model(env::Env, name::String, sense::Symbol)
    @assert isascii(name)
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

Base.unsafe_convert(ty::Type{Ptr{Cvoid}}, model::Model) = model.ptr_model::Ptr{Cvoid}

function free_model(model::Model)
    if model.ptr_model != C_NULL
        @grb_ccall(freemodel, Nothing, (Ptr{Cvoid},), model.ptr_model)
        model.ptr_model = C_NULL
    end
end

function copy(model::Model)
    pm::Ptr{Cvoid} = C_NULL
    if model.ptr_model != C_NULL
        pm = @grb_ccall(copymodel, Ptr{Cvoid}, (Ptr{Cvoid},), model.ptr_model)
        if pm == C_NULL
            error("Failed to copy a Gurobi model.")
        end
    end
    Model(model.env, pm)
end

function update_model!(model::Model)
    @assert model.ptr_model != C_NULL
    ret = @grb_ccall(updatemodel, Cint, (Ptr{Cvoid},), model.ptr_model)
    if ret != 0
        throw(GurobiError(model.env, ret))
    end
    nothing
end

function reset_model!(model::Model)
    @assert model.ptr_model != C_NULL
    ret = @grb_ccall(resetmodel, Cint, (Ptr{Cvoid},), model.ptr_model)
    if ret != 0
        throw(GurobiError(model.env, ret))
    end
    nothing
end

# read / write file

function read(model::Model, filename::String)
    @assert isascii(filename) # TODO: support non-ascii file names
        ret = @grb_ccall(read, Cint,
        (Ptr{Cvoid}, Ptr{UInt8}),
        model.ptr_model, filename)
    if ret != 0
        throw(GurobiError(model.env, ret))
    end
    nothing
end

function read_model(model::Model, filename::String)
    @assert isascii(filename) # TODO: support non-ascii file names
    @assert is_valid(model.env)
    a = Ref{Ptr{Cvoid}}()
    ret = @grb_ccall(readmodel, Cint,
        (Ptr{Cvoid}, Ptr{UInt8}, Ptr{Ptr{Cvoid}}),
        model.env, filename, a)
    if ret != 0
        throw(GurobiError(model.env, ret))
    end
    model.ptr_model = a[]
    nothing
end

function write_model(model::Model, filename::String)
    @assert isascii(filename) # TODO: support non-ascii file names
    ret = @grb_ccall(write, Cint, (Ptr{Cvoid}, Ptr{UInt8}),
        model.ptr_model, filename)
    if ret != 0
        throw(GurobiError(model.env, ret))
    end
    nothing
end

function tune_model(model::Model)
    ret = @grb_ccall(tunemodel, Cint, (Ptr{Cvoid},), model.ptr_model)
    if ret != 0
        throw(GurobiError(model.env, ret))
    end
    nothing
end

function get_tune_result!(model::Model,i::Int)
    ret = @grb_ccall(gettuneresult, Cint, (Ptr{Cvoid}, Cint), model.ptr_model, i)
    if ret != 0
        throw(GurobiError(model.env, ret))
    end
    nothing
end

terminate(model::Model) = @grb_ccall(terminate, Nothing, (Ptr{Cvoid},), model.ptr_model)

# Presolve the model but don't solve. For some reason this is not
# documented for the C interface, but is for all the other interfaces.
# Source: https://twitter.com/iaindunning/status/519620465992556544
function presolve_model(model::Model)
    ret = @grb_ccall(presolvemodel, Ptr{Cvoid}, (Ptr{Cvoid},), model.ptr_model)
    if ret == C_NULL
        # Presumably failed to return a model
        error("presolve_model did not return a model")
    end
    return Model(model.env, ret)
end

function fixed_model(model::Model)
    @assert model.ptr_model != C_NULL
    fixed::Ptr{Cvoid} = C_NULL
    fixed = @grb_ccall(fixedmodel, Ptr{Cvoid}, (Ptr{Cvoid},), model.ptr_model)
    if fixed == C_NULL
        error("Unable to create fixed model")
    end
    Model(model.env, fixed)
end
