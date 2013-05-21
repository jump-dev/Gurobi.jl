# GurobiSolverInterface
# Standardized LP interface

export GurobiSolver,
    model,
    loadproblem,
    writeproblem,
    getvarLB,
    setvarLB,
    getvarLB,
    setvarLB,
    getconstrLB,
    setconstrLB,
    getconstrUB,
    setconstrUB,
    getobj,
    setobj,
    addvar,
    addconstr,
    updatemodel,
    setsense,
    getsense,
    numvar,
    numconstr,
    optimize,
    status,
    getobjval,
    getobjbound,
    getsolution,
    getconstrsolution,
    getreducedcosts,
    getconstrduals,
    getrawsolver,
    setvartype

type GurobiSolver <: LinprogSolver
	inner::Model
end

# Not right...
function model(;options...)
	if length(options) != 0; warn("Options not yet supported"); end
	m = GurobiSolver(gurobi_model(Env(),""))
	return m
end

# Creates new env? doesn't seem right
function loadproblem(m::GurobiSolver, filename::String)
  m.inner = read_model(Env(), filename)
end

function loadproblem(m::GurobiSolver, A, collb, colub, obj, rowlb, rowub)
  reset_model!(m.inner)
  add_cvars!(m.inner, float(obj), float(collb), float(colub))
  update_model!(m.inner)
  add_rangeconstrs!(m.inner, A, float(rowlb), float(rowub))
  update_model!(m.inner)
end

function writeproblem(m::GurobiSolver, filename::String)
  write_model(m.inner, filename)
end


######
# More functions need to be exposed here really
######

# TODO
function updatemodel(m::GurobiSolver)
  error("Not Implemented - what is this?")
end

function setsense(m::GurobiSolver,sense)
  if sense == :Min
    set_sense!(m.inner, :minimize)
  elseif sense == :Max
    set_sense!(m.inner, :maximize)
  else
    error("Unrecognized objective sense $sense")
  end
end
function getsense(m::GurobiSolver)
  v = get_int_attr(m.inner, "ModelSense")
  if v == -1 
    return :Max 
  else
    return :Min
  end
end

numvar(m::GurobiSolver)    = num_vars(m.inner)
numconstr(m::GurobiSolver) = num_constrs(m.inner)

optimize(m::GurobiSolver)  = optimize(m.inner)

function status(m::GurobiSolver)
  s = get_status(m.inner)
  if s == :optimal
    return :Optimal
  elseif s == :infeasible
    return :Infeasible
  elseif s == :inf_or_unbd
    return :Unbounded
  elseif s == :iteration_limit || s == :node_limit || s == :time_limit || s == :solution_limit
    return :UserLimit
  else
    error("Internal library error")
  end
end

getobjval(m::GurobiSolver)   = get_objval(m.inner)
getobjbound(m::GurobiSolver) = get_objbound(m.inner)
getsolution(m::GurobiSolver) = get_solution(m.inner)

# TODO
function getconstrsolution(m::GurobiSolver)
  error("GurobiSolver: Not implemented (need to do Ax manually?)")
end

getreducedcosts(m::GurobiSolver) = get_dbl_attrarray(m.inner, "RC", 1, num_vars(m.inner))
getconstrduals(m::GurobiSolver)  = get_dbl_attrarray(m.inner, "Pi", 1, num_constrs(m.inner))

getrawsolver(m::GurobiSolver) = m.inner

setvartype(m::GurobiSolver, vartype) =
    set_char_attr_array!(m.inner, "VType", 1, length(vartype), vartype)
