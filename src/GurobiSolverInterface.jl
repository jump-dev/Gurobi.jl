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
    getsolution,
    getconstrsolution,
    getreducedcosts,
    getconstrduals,
    getrawsolver

type GurobiSolver <: LinprogSolver
	inner::Model
end

# Not right...
function model()
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

function optimize(m::GurobiSolver)
  optimize(m.inner)
end

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

# TODO
function getreducedcosts(m::GurobiSolver)
  return zeros(num_vars(model))
end

# TODO
function getconstrduals(m::GurobiSolver)
  return zeros(num_constrs(model))
end

function getrawsolver(m::GurobiSolver)
  return m.inner
end
