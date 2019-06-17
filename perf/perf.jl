using Gurobi

function print_help()
    println("""
    Usage
        perf.jl [arg] [name]

    [arg]
        --new       Begin a new benchmark comparison
        --compare   Run another benchmark and compare to existing

    [name]          A name for the benchmark test

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
    const Benchmarks = Gurobi.MOI.Benchmarks
    const GUROBI_ENV = Gurobi.Env()
    const suite = Benchmarks.suite() do
        Gurobi.Optimizer(GUROBI_ENV, OutputFlag = 0)
    end
    if ARGS[1] == "--new"
        Benchmarks.create_baseline(
            suite, ARGS[2]; directory = @__DIR__, verbose = true
        )
    elseif ARGS[1] == "--compare"
        Benchmarks.compare_against_baseline(
            suite, ARGS[2]; directory = @__DIR__, verbose = true
        )
    else
        print_help()
    end
end
