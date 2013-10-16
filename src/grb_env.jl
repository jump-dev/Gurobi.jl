# Gurobi environment and other supporting facilities


type Env
    ptr_env::Ptr{Void}
    
    function Env()
        a = Array(Ptr{Void}, 1)
        ret = @grb_ccall(loadenv, Cint, (Ptr{Ptr{Void}}, Ptr{Uint8}), 
            a, C_NULL)
        if ret != 0
            error("Failed to create environment (error code = $ret).")
        end
        env = new(a[1])
        # finalizer(env, free_env)  ## temporary disable: which tends to sometimes caused warnings
        env
    end
end

convert(ty::Type{Ptr{Void}}, env::Env) = env.ptr_env::Ptr{Void}

function is_valid(env::Env)
    env.ptr_env != C_NULL
end

function free_env(env::Env)
    if env.ptr_env != C_NULL
        @grb_ccall(freeenv, Void, (Ptr{Void},), env.ptr_env)
        env.ptr_env = C_NULL
    end
end

function get_error_msg(env::Env)
    @assert env.ptr_env != C_NULL
    sz = @grb_ccall(geterrormsg, Ptr{Uint8}, (Ptr{Void},), env.ptr_env)
    bytestring(sz)
end

# error

type GurobiError
    code::Int
    msg::ASCIIString 
    
    function GurobiError(env::Env, code::Integer)
        new(convert(Int, code), get_error_msg(env))
    end
end

