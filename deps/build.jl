using BinDeps

@BinDeps.setup

libgurobi = library_dependency("libgurobi",aliases=["libgurobi51","libgurobi55"])

provides(Binaries, joinpath(ENV["GUROBI_HOME"],"lib"), libgurobi)
