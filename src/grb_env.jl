# Gurobi environment and other supporting facilities


type Env
    ptr_env::Ptr{Void}

    """
        Env(silent::Bool=false)

    Construct a Gurobi environment. The environment maintains parameters of the
    Gurobi solver and its license. In particular, constructing an environment
    will obtain a new Gurobi license token. If you are operating in a setting
    where the number of concurrent Gurobi tokens is limited, you may want to
    construct a single environment ahead of time and share it among all the
    models solved by your program.

    The `silent` option temporarily redirects the standard output stream
    while constructing the environment. This may be useful in preventing
    Gurobi from printing messages to the screen every time an environment
    is constructed.
    """
    function Env(silent::Bool=false)
        a = Array{Ptr{Void}}(1)
        if silent
            _stdout = STDOUT
            stdout_rd, stdout_wr = redirect_stdout()
            ret = try
                @grb_ccall(loadenv, Cint, (Ptr{Ptr{Void}}, Ptr{UInt8}),
                    a, C_NULL)
            finally
                redirect_stdout(_stdout)
                close(stdout_rd)
                close(stdout_wr)
            end
        else
            ret = @grb_ccall(loadenv, Cint, (Ptr{Ptr{Void}}, Ptr{UInt8}),
                a, C_NULL)
        end

        if ret != 0
            if ret == 10009
                error("Invalid Gurobi license")
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
    unsafe_string(sz)
end

# error

type GurobiError
    code::Int
    msg::String

    function GurobiError(env::Env, code::Integer)
        new(convert(Int, code), get_error_msg(env))
    end
end

