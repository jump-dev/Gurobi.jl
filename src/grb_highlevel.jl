# High level model construction

function gurobi_model(env::Env;
	name::ASCIIString="", 
	f::Vector{Float64}=Float64[],
    sense::Symbol=:minimize,       # :minimize or :maximize
    H=nothing,                     # quadratic matrix 
    A::Union(ConstrMat, Nothing)=nothing,         # inequality constraints
    b::Union(Vector{Float64}, Nothing)=nothing,           
    Aeq::Union(ConstrMat, Nothing)=nothing,       # equality constraints 
    beq::Union(Vector{Float64}, Nothing)=nothing, 
    lb::Bounds=nothing,    # upper bounds 
    ub::Bounds=nothing)    # lower bounds

	# check f
	if isempty(f)
		error("f must be specified.")
	end

    # create model
    model = Model(env, name)
    
    # set sense
    if sense != :minimize
        set_sense!(model, sense)
    end

    # add variables
    add_cvars!(model, f, lb, ub)
    update_model!(model)
    
    # add qpterms
    if !is(H, nothing)
        add_qpterms!(model, H)
    end
    
    # add constraints
    if !is(A, nothing) && !is(b, nothing)
        add_constrs!(model, A, '<', b)
    end
    
    if !is(Aeq, nothing) && !is(beq, nothing)
        add_constrs!(model, Aeq, '=', beq)
    end
    update_model!(model)
    
    return model
end


