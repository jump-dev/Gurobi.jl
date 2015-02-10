# Gurobi callbacks

type CallbackData
    cbdata::Ptr{Void}
    model::Model
end

function gurobi_callback_wrapper(ptr_model::Ptr{Void}, cbdata::Ptr{Void}, where::Cint, userdata::Ptr{Void})

    callback,model = unsafe_pointer_to_objref(userdata)::(Function,Model)
    callback(CallbackData(cbdata,model), where)
    return convert(Cint,0)
end

# User callback function should be of the form:
# callback(cbdata::CallbackData, where::Cint)

function set_callback_func!(model::Model, callback::Function)
    
    grbcallback = cfunction(gurobi_callback_wrapper, Cint, (Ptr{Void}, Ptr{Void}, Cint, Ptr{Void}))
    usrdata = (callback,model)
    ret = @grb_ccall(setcallbackfunc, Cint, (Ptr{Void}, Ptr{Void}, Any), model.ptr_model, grbcallback, usrdata)
    if ret != 0
        throw(GurobiError(model.env, ret))
    end
    # we need to keep a reference to the callback function
    # so that it isn't garbage collected
    model.callback = usrdata
    nothing
end

export CallbackData, set_callback_func!

for f in (:cbcut, :cblazy)
    @eval function ($f)(cbdata::CallbackData, ind::Vector{Cint}, val::Vector{Float64},
        sense::Char, rhs::Float64)

        len = length(ind)
        @assert length(val) == len

        ret = @grb_ccall($f, Cint, (Ptr{Void},Cint,Ptr{Cint},Ptr{Float64},
            Char,Float64), cbdata.cbdata, len, ind.-1, val, sense, rhs)
        if ret != 0
            throw(GurobiError(cbdata.model.env, ret))
        end
    end
end

export cbcut, cblazy

function cbsolution(cbdata::CallbackData, sol::Vector{Float64})
    nvar = num_vars(cbdata.model)
    @assert length(sol) >= nvar

    ret = @grb_ccall(cbsolution, Cint, (Ptr{Void},Ptr{Float64}),
        cbdata.cbdata, sol)
    if ret != 0
        throw(GurobiError(cbdata.model.env, ret))
    end
end


function cbget{T}(::Type{T},cbdata::CallbackData, where::Cint, what::Integer)

    out = Array(T,1)
    ret = @grb_ccall(cbget, Cint, (Ptr{Void}, Cint, Cint, Ptr{T}),
        cbdata.cbdata, where, convert(Cint,what), out)
    if ret != 0
        throw(GurobiError(cbdata.model.env, ret))
    end
    return out[1]
end


# Callback constants
# grep GRB_CB gurobi_c.h | awk '{ print "const " substr($2,5) " = " $3; }'
const CB_POLLING = 0
const CB_PRESOLVE = 1
const CB_SIMPLEX = 2
const CB_MIP = 3
const CB_MIPSOL = 4
const CB_MIPNODE = 5
const CB_MESSAGE = 6
const CB_BARRIER = 7

export CB_POLLING, CB_PRESOLVE, CB_SIMPLEX, CB_MIP,
       CB_MIPSOL, CB_MIPNODE, CB_MESSAGE, CB_BARRIER

# grep GRB_CB gurobi_c.h | awk '{ print "(\"" tolower(substr($2,8)) "\"," $3 ")"; }'
const cbconstants = [
("pre_coldel",1000,Cint),
("pre_rowdel",1001,Cint),
("pre_senchg",1002,Cint),
("pre_bndchg",1003,Cint),
("pre_coechg",1004,Cint),
("spx_itrcnt",2000,Float64),
("spx_objval",2001,Float64),
("spx_priminf",2002,Float64),
("spx_dualinf",2003,Float64),
("spx_ispert",2004,Float64),
("mip_objbst",3000,Float64),
("mip_objbnd",3001,Float64),
("mip_nodcnt",3002,Float64),
("mip_solcnt",3003,Cint),
("mip_cutcnt",3004,Cint),
("mip_nodlft",3005,Float64),
("mip_itrcnt",3006,Float64),
###("mipsol_sol",4001),
("mipsol_obj",4002,Float64),
("mipsol_objbst",4003,Float64),
("mipsol_objbnd",4004,Float64),
("mipsol_nodcnt",4005,Float64),
("mipsol_solcnt",4006,Cint),
("mipnode_status",5001,Cint),
###("mipnode_rel",5002),
("mipnode_objbst",5003,Float64),
("mipnode_objbnd",5004,Float64),
("mipnode_nodcnt",5005,Float64),
("mipnode_solcnt",5006,Cint),
##("mipnode_brvar",5007), -- undocumented
##("msg_string",6001), -- not yet implemented: 
### documentation is unclear on output type
("runtime",6002, Float64),
("barrier_itrcnt",7001,Cint),
("barrier_primobj",7002,Float64),
("barrier_dualobj",7003,Float64),
("barrier_priminf",7004,Float64),
("barrier_dualinf",7005,Float64),
("barrier_compl",7006,Float64)]

for (cname,what,T) in cbconstants
    fname = symbol("cbget_$cname")
    @eval ($fname)(cbdata::CallbackData, where::Cint) = cbget($T, cbdata, where, $what)
    eval(Expr(:export,fname))
end

for (fname, what) in ((:cbget_mipsol_sol, 4001), (:cbget_mipnode_rel, 5002))
    @eval function ($fname)(cbdata::CallbackData, where::Cint, out::Vector{Float64})
        nvar = num_vars(cbdata.model)
        @assert length(out) >= nvar
        ret = @grb_ccall(cbget, Cint, (Ptr{Void}, Cint, Cint, Ptr{Float64}),
                         cbdata.cbdata, where, $what, out)
        if ret != 0
            throw(GurobiError(cbdata.model.env, ret))
        end
    end
    @eval function ($fname)(cbdata::CallbackData, where::Cint)
        nvar = num_vars(cbdata.model)
        out = Array(Float64,nvar)
        ($fname)(cbdata, where, out)
        return out
    end
    eval(Expr(:export,fname))
end



