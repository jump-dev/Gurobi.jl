# Gurobi environment and other supporting facilities


type Env
    ptr_env::Ptr{Void}

    function Env(cloudhost::ASCIIString="", cloudpassword::ASCIIString="")
        a = Array(Ptr{Void}, 1)
        if length(cloudhost) > 0 # cloud server
            ret = @grb_ccall(loadclientenv, Cint,
                (Ptr{Ptr{Void}}, Ptr{UInt8}, Ptr{UInt8}, Cint, Ptr{UInt8}, Cint, Cdouble),
                a, C_NULL, cloudhost, -1, cloudpassword, 0, -1)
        else # local installation
            ret = @grb_ccall(loadenv, Cint, (Ptr{Ptr{Void}}, Ptr{UInt8}),
                a, C_NULL)
        end
        if ret != 0
            if ret == 10009
                error("Invalid Gurobi license")
            elseif ret == 10022
                error("Problem communicating with the Gurobi Compute Server")
            else
                error("Failed to create environment (error $ret).")
            end
        end
        env = new(a[1])
        # finalizer(env, free_env)  ## temporary disable: which tends to sometimes caused warnings
        env
    end
end

Base.unsafe_convert(ty::Type{Ptr{Void}}, env::Env) = env.ptr_env::Ptr{Void}

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
    sz = @grb_ccall(geterrormsg, Ptr{UInt8}, (Ptr{Void},), env.ptr_env)
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

