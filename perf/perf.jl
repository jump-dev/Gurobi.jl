using Gurobi, BenchmarkTools

const MOI = Gurobi.MOI
const GUROBI_ENV = Gurobi.Env()
function new_model()
    return Gurobi.Optimizer(GUROBI_ENV, OutputFlag=0)
end

suite = BenchmarkGroup()

function add_variable()
    model = new_model()
    for i in 1:10_000
        MOI.add_variable(model)
    end
    return model
end
suite["add_variable"] = @benchmarkable add_variable()

function add_variables()
    model = new_model()
    MOI.add_variables(model, 10_000)
    return model
end
suite["add_variables"] = @benchmarkable add_variables()

function add_variable_constraint()
    model = new_model()
    x = MOI.add_variables(model, 10_000)
    for (i, xi) in enumerate(x)
        MOI.add_constraint(model, MOI.SingleVariable(xi), MOI.LessThan(1.0 * i))
    end
    return model
end
suite["add_variable_constraint"] = @benchmarkable add_variable_constraint()

function add_variable_constraints()
    model = new_model()
    x = MOI.add_variables(model, 10_000)
    MOI.add_constraints(
        model,
        MOI.SingleVariable.(x),
        MOI.LessThan.(1.0:10_000.0)
    )
    return model
end
suite["add_variable_constraints"] = @benchmarkable add_variable_constraints()

function delete_variable()
    model = new_model()
    x = MOI.add_variables(model, 1_000)
    MOI.add_constraint.(model, MOI.SingleVariable.(x), Ref(MOI.LessThan(1.0)))
    MOI.delete.(model, x)
    return model
end
suite["delete_variable"] = @benchmarkable delete_variable()

function delete_variable_constraint()
    model = new_model()
    x = MOI.add_variables(model, 1_000)
    cons = MOI.add_constraint.(model, MOI.SingleVariable.(x), Ref(MOI.LessThan(1.0)))
    for con in cons
        MOI.delete(model, con)
    end
    cons = MOI.add_constraint.(model, MOI.SingleVariable.(x), Ref(MOI.LessThan(1.0)))
    MOI.set(model, MOI.ObjectiveSense(), MOI.MAX_SENSE)
    MOI.set(model,
        MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(),
        MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(1.0, x), 0.0)
    )
    MOI.optimize!(model)
    for con in cons
        MOI.delete(model, con)
    end
    return model
end
suite["delete_variable_constraint"] = @benchmarkable delete_variable_constraint()

function add_constraint()
    model = new_model()
    index = MOI.add_variables(model, 10_000)
    for (i, x) in enumerate(index)
        MOI.add_constraint(
            model,
            MOI.ScalarAffineFunction([MOI.ScalarAffineTerm(1.0, x)], 0.0),
            MOI.LessThan(1.0 * i)
        )
    end
    return model
end
suite["add_constraint"] = @benchmarkable add_constraint()

function add_constraints()
    model = new_model()
    x = MOI.add_variables(model, 10_000)
    MOI.add_constraints(
        model,
        [MOI.ScalarAffineFunction([MOI.ScalarAffineTerm(1.0, xi)], 0.0) for xi in x],
        MOI.LessThan.(1:1.0:10_000)
    )
    return model
end
suite["add_constraints"] = @benchmarkable add_constraints()

function delete_constraint()
    model = new_model()
    index = MOI.add_variables(model, 1_000)
    cons = Vector{
        MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64}, MOI.LessThan{Float64}}
    }(undef, 1_000)
    for (i, x) in enumerate(index)
        cons[i] = MOI.add_constraint(
            model,
            MOI.ScalarAffineFunction([MOI.ScalarAffineTerm(1.0, x)], 0.0),
            MOI.LessThan(1.0 * i)
        )
    end
    for con in cons
        MOI.delete(model, con)
    end
    return model
end
suite["delete_constraint"] = @benchmarkable delete_constraint()

function new_test(suite, params_filename, results_filename)
    tune!(suite)
    BenchmarkTools.save(params_filename, params(suite))
    results = run(suite, verbose = true)
    BenchmarkTools.save(results_filename, results)
end

function compare_test(suite, params_filename, results_filename)
    if !isfile(params_filename) || !isfile(results_filename)
        error("You must run with --new first.")
    end
    loadparams!(
        suite,
        BenchmarkTools.load(params_filename)[1],
        :evals, :samples
    )
    new_results = run(suite, verbose = true)
    old_results = BenchmarkTools.load(results_filename)[1]
    println(stdout)
    println("========== Results ==========")
    println(stdout)
    for key in keys(new_results)
        judgement = judge(
            BenchmarkTools.median(new_results[key]),
            BenchmarkTools.median(old_results[key])
        )
        println(stdout, key)
        show(stdout, MIME"text/plain"(), judgement)
        println(stdout)
    end
end

function print_help()
    println("""
    Usage
        perf.jl [arg] [name]

    [arg]
        --new       Begin a new benchmark comparison
        --compare   Run another benchmark and compare to existing

    [name]          A name for the benchmark test. Will create a folder in /perf.

    Examples
        git checkout master
        julia perf.jl --new master
        git checkout approach_1
        julia perf.jl --new approach_1
        git checkout approach_2
        julia perf.jl --compare master
        julia perf.jl --compare approach_1
    """)
end

if length(ARGS) != 2
    print_help()
else
    test_path = joinpath(@__DIR__, ARGS[2])
    if !isdir(test_path)
        mkdir(test_path)
    end
    params_filename = joinpath(test_path, "params.json")
    results_filename = joinpath(test_path, "results.json")
    if ARGS[1] == "--new"
        new_test(suite, params_filename, results_filename)
    elseif ARGS[1] == "--compare"
        compare_test(suite, params_filename, results_filename)
    else
        print_help()
    end
end
