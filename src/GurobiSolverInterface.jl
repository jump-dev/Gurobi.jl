# GurobiMathProgInterface
# Standardized LP interface

export GurobiSolver

type GurobiMathProgModel <: AbstractMathProgModel
	inner::Model
end

immutable GurobiSolver <: AbstractMathProgSolver
    options
end
GurobiSolver(;kwargs...) = GurobiSolver(kwargs)

function GurobiMathProgModel(;options...)
	env = Env()
	for (name,value) in options
		setparam!(env, string(name), value)
	end
	m = GurobiMathProgModel(Model(env,""))
	return m
end

model(s::GurobiSolver) = GurobiMathProgModel(;s.options...)

# Creates new env? doesn't seem right
#function loadproblem!(m::GurobiMathProgModel, filename::String)
#  m.inner = read_model(Env(), filename)
#end

function loadproblem!(m::GurobiMathProgModel, A, collb, colub, obj, rowlb, rowub, sense)
  reset_model!(m.inner)
  add_cvars!(m.inner, float(obj), float(collb), float(colub))
  update_model!(m.inner)
  add_rangeconstrs!(m.inner, float(A), float(rowlb), float(rowub))
  update_model!(m.inner)
  setsense!(m,sense)
end

function writeproblem(m::GurobiMathProgModel, filename::String)
  write_model(m.inner, filename)
end


######
# Incomplete implementation of AbstractMathProgModel
######

function updatemodel!(m::GurobiMathProgModel)
  update_model!(m.inner)
end

function setsense!(m::GurobiMathProgModel,sense)
  if sense == :Min
    set_sense!(m.inner, :minimize)
  elseif sense == :Max
    set_sense!(m.inner, :maximize)
  else
    error("Unrecognized objective sense $sense")
  end
end
function getsense(m::GurobiMathProgModel)
  v = get_intattr(m.inner, "ModelSense")
  if v == -1 
    return :Max 
  else
    return :Min
  end
end

numvar(m::GurobiMathProgModel)    = num_vars(m.inner)
numconstr(m::GurobiMathProgModel) = num_constrs(m.inner)

optimize!(m::GurobiMathProgModel)  = optimize(m.inner)

function status(m::GurobiMathProgModel)
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

getobjval(m::GurobiMathProgModel)   = get_objval(m.inner)
getobjbound(m::GurobiMathProgModel) = get_objbound(m.inner)
getsolution(m::GurobiMathProgModel) = get_solution(m.inner)

# TODO
function getconstrsolution(m::GurobiMathProgModel)
  error("GurobiMathProgModel: Not implemented (need to do Ax manually?)")
end

getreducedcosts(m::GurobiMathProgModel) = get_dblattrarray(m.inner, "RC", 1, num_vars(m.inner))
getconstrduals(m::GurobiMathProgModel)  = get_dblattrarray(m.inner, "Pi", 1, num_constrs(m.inner))

getrawsolver(m::GurobiMathProgModel) = m.inner

setvartype!(m::GurobiMathProgModel, vartype) =
    set_charattrarray!(m.inner, "VType", 1, length(vartype), vartype)
