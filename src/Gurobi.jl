module Gurobi
   
using BinDeps
@BinDeps.load_dependencies

    ### imports

    import Base.convert, Base.show, Base.copy

    # Standard LP interface
    require(joinpath(Pkg.dir("MathProgBase"),"src","LinprogSolverInterface.jl"))
    importall LinprogSolverInterface

    ### exports
    export 

    # grb_env
    free_env,

    # grb_params
    getparam, setparam!, setparams!,

    # grb_model
    set_sense!, update_model!, reset_model!, 
    read_model, write_model,

    # grb_vars
    GRB_CONTINUOUS, GRB_BINARY, GRB_INTEGER,
    add_var!, add_vars!, add_cvar!, add_cvars!,
    add_bvar!, add_bvars!, add_ivar!, add_ivars!,

    # grb_constrs
    add_constr!, add_constrs!, add_rangeconstr!, add_rangeconstrs!,

    # grb_quad
    add_qpterms!, add_qconstr!,

    # higher level
    gurobi_model, qp_model,

    # grb_solve
    optimize, computeIIS, get_solution,
    get_status, OptimInfo, get_optim_info, get_objval


    ### common support
    
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


    ### include source files

    include("grb_env.jl")
    include("grb_params.jl")

    include("grb_model.jl")
    include("grb_vars.jl")
    include("grb_attrs.jl")
    include("grb_constrs.jl")
    include("grb_quad.jl")
    include("grb_highlevel.jl")

    include("grb_solve.jl")
    include("grb_callbacks.jl")

    include("GurobiSolverInterface.jl")

    include("deprecates.jl")
end
