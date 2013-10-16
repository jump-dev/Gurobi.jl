module Gurobi
   
using BinDeps
@BinDeps.load_dependencies

    # Standard LP interface
    require(joinpath(Pkg.dir("MathProgBase"),"src","LinprogSolverInterface.jl"))
    importall LinprogSolverInterface

    export getparam, setparam, setparams!

    export update_model!, reset_model!
    export read_model, write_model
    export set_sense!
    export add_var!, add_vars!, add_cvar!, add_cvars!
    export add_bvar!, add_bvars!, add_ivar!, add_ivars!
    export add_constr!, add_constrs!, add_qpterms!, add_qconstr!
    export add_rangeconstr!, add_rangeconstrs!
    
    export qp_model
    export optimize, computeIIS
    
    export get_status, OptimInfo, get_optim_info, get_objval
    export get_solution
    
    import Base.convert, Base.show, Base.copy
    
    # macro to call a Gurobi C function
    macro grb_ccall(func, args...)
        f = "GRB$(func)"
        quote
            ccall(($f,libgurobi), $(args...))
        end
    end

    # Gurobi library version
    function getlibversion()
        _major = Cint[0]
        _minor = Cint[0]
        _tech = Cint[0]
        @grb_ccall(version, Void, (Ptr{Cint}, Ptr{Cint}, Ptr{Cint}), _major, _minor, _tech)
        return VersionNumber(_major[1], _minor[1], _tech[1])        
    end

    # version need not be export
    # one can write Gurobi.version to get the version numbers
    const version = getlibversion()

    # include source files

    include("grb_env.jl")
    include("grb_params.jl")
    include("grb_model.jl")
    include("grb_solve.jl")
    include("grb_callbacks.jl")

    include("GurobiSolverInterface.jl")

    include("deprecates.jl")
end
