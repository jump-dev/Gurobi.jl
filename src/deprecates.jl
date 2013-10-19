
# TODO : to be deprecated in favor of Gurobi.Model
function qp_model(env::Env, name::ASCIIString, 
    H::Union(Vector{Float64}, Matrix{Float64}, SparseMatrixCSC{Float64}, Float64), 
    f::Vector{Float64}, 
    A::Union(CoeffMat, Nothing), 
    b::Union(Vector{Float64}, Nothing), 
    Aeq::Union(CoeffMat, Nothing), 
    beq::Union(Vector{Float64}, Nothing), 
    lb::Bounds, ub::Bounds)

    Base.warn_once("qp_model is to be deprecated in favor of gurobi_model with keyword arguments.")
    
    # create model
    model = Model(env, name)
    
    # add variables
    add_cvars!(model, f, lb, ub)
    update_model!(model)
    
    # add qpterms
    
    _add_qpterms!(model, H)
    
    # add constraints
    if A != nothing && b != nothing
        add_constrs!(model, A, '<', b)
    end
    
    if Aeq != nothing && beq != nothing
        add_constrs!(model, Aeq, '=', beq)
    end
    update_model!(model)
    
    model
end

qp_model(env::Env, name, H, f, A, b, Aeq, beq) = qp_model(env, name, H, f, A, b, Aeq, beq, nothing, nothing)
qp_model(env::Env, name, H, f, A, b) = qp_model(env, name, H, f, A, b, nothing, nothing, nothing, nothing)
