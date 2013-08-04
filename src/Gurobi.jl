module Gurobi
   
using BinDeps
@BinDeps.load_dependencies

	# Standard LP interface
	require(joinpath(Pkg.dir("MathProgBase"),"src","LinprogSolverInterface.jl"))
    importall LinprogSolverInterface

    export set_int_param!, set_dbl_param!, get_int_param, get_dbl_param

    export gurobi_model, update_model!, reset_model!
    export read_model, write_model
    export set_sense!
    export add_var!, add_vars!, add_cvar!, add_cvars!
    export add_bvar!, add_bvars!, add_ivar!, add_ivars!
    export add_constr!, add_constrs!, add_qpterms!, add_qconstr!
	export add_rangeconstr!, add_rangeconstrs!
    
    export lp_model, qp_model
    export optimize, computeIIS
    
    export get_status, OptimInfo, get_optim_info, get_objval
    export get_solution
    
    import Base.convert, Base.show, Base.copy
    
    macro grb_ccall(func, args...)
        f = "GRB$(func)"
        quote
            ccall(($f,libgurobi), $(args...))
        end
    end

    include("grb_env.jl")
    include("grb_model.jl")
    include("grb_solve.jl")

	include("GurobiSolverInterface.jl")
end
