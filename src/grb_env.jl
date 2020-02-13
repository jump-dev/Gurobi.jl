# Gurobi environment and other supporting facilities

mutable struct Env
    ptr_env::Ptr{Cvoid}
end

# Note: `Env()` is an outer constructor because sometimes we want to
# make an `Env` object from a pointer returned from `GRBgetenv`, and
# sometimes we want to create a new `Env`.
function Env()
    a = Ref{Ptr{Cvoid}}()
    ret = @grb_ccall(
        loadenv,
        Cint,
        (Ptr{Ptr{Cvoid}}, Ptr{UInt8}),
        a, C_NULL
    )
    if ret == 10009
        error("Invalid Gurobi license")
    elseif ret != 0
        error("Failed to create environment (error $ret).")
    end
    # TODO(odow): no finalizer is set for Env, because when a model + Env
    # falls out of scope, Julia's GC will sometimes GC the env first, and
    # then the model. This causes an error. Thus, `Env`'s are finalized in
    # the `Model` constructor.
    # We should probably add a flag to differentiate between envs created
    # manually, and envs created by `Model()`. CPLEX.jl does something
    # similar.
    # finalizer(env, free_env)
    return Env(a[])
end

Base.unsafe_convert(ty::Type{Ptr{Cvoid}}, env::Env) = env.ptr_env::Ptr{Cvoid}

function is_valid(env::Env)
    env.ptr_env != C_NULL
end

function free_env(env::Env)
    if env.ptr_env != C_NULL
        @grb_ccall(freeenv, Nothing, (Ptr{Cvoid},), env.ptr_env)
        env.ptr_env = C_NULL
    end
end

function get_error_msg(env::Env)
    @assert env.ptr_env != C_NULL
    sz = @grb_ccall(geterrormsg, Ptr{UInt8}, (Ptr{Cvoid},), env.ptr_env)
    unsafe_string(sz)
end

# error

mutable struct GurobiError <: Exception
    code::Int
    msg::String

    function GurobiError(env::Env, code::Integer)
        new(convert(Int, code), get_error_msg(env))
    end
end
