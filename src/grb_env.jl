# Gurobi environment and other supporting facilities


mutable struct Env
    ptr_env::Ptr{Cvoid}

    function Env()
        a = Ref{Ptr{Cvoid}}()
        ret = @grb_ccall(loadenv, Cint, (Ptr{Ptr{Cvoid}}, Ptr{UInt8}),
            a, C_NULL)
        if ret != 0
            if ret == 10009
                error("Invalid Gurobi license")
            else
                error("Failed to create environment (error $ret).")
            end
        end
        env = new(a[])
        # finalizer(env, free_env)  ## temporary disable: which tends to sometimes caused warnings
        env
    end
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
