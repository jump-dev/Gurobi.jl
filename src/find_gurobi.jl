# Locating Gurobi library
#
# Current approach simply looks at an environment variable GUROBI_LIB
#

function _find_library_path()
	if !has(ENV, "GUROBI_LIB")
		error("Please set the environment variable GUROBI_LIB to specify the dynamic library path.")
	end
	p = ENV["GUROBI_LIB"]
    if !isfile(p)
        error("Cannot find library file $p.")
    end
    p
end
	
const libgurobi_path = _find_library_path()
const lib_gurobi = dlopen(libgurobi_path, RTLD_GLOBAL)

# function pointers

macro grb_def(fun)
    local fptr = symbol(string(fun, "_ptr"))
    @eval begin
        $(fptr) = C_NULL
        function $(fun)()
            global $(fptr)
            if $(fptr) == C_NULL
                $(fptr) = dlsym(lib_gurobi, $(Meta.quot(fun)))
            end
            $(fptr)::Ptr{Void}
        end
    end
end

@grb_def GRBloadenv
@grb_def GRBfreeenv
@grb_def GRBgeterrormsg

@grb_def GRBnewmodel
@grb_def GRBcopymodel
@grb_def GRBfreemodel
@grb_def GRBupdatemodel
@grb_def GRBresetmodel

@grb_def GRBreadmodel
@grb_def GRBwrite

@grb_def GRBaddvar
@grb_def GRBaddvars
@grb_def GRBaddconstr
@grb_def GRBaddconstrs
@grb_def GRBaddqpterms

@grb_def GRBgetintattr
@grb_def GRBgetdblattr
@grb_def GRBgetstrattr
@grb_def GRBgetintattrarray
@grb_def GRBgetdblattrarray
@grb_def GRBgetstrattrarray

@grb_def GRBsetintattr
@grb_def GRBsetdblattr
@grb_def GRBsetstrattr

@grb_def GRBgetintparam
@grb_def GRBgetdblparam
@grb_def GRBgetstrparam
@grb_def GRBsetintparam
@grb_def GRBsetdblparam
@grb_def GRBsetstrparam

@grb_def GRBoptimize
