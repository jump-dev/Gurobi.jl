using BinDeps

@BinDeps.setup

libgurobi = library_dependency("libgurobi",aliases=["libgurobi51","gurobi51","libgurobi55","gurobi55","libgurobi56","gurobi56"])

@unix_only provides(Binaries, joinpath(ENV["GUROBI_HOME"],"lib"), libgurobi)
@windows_only provides(Binaries, joinpath(ENV["GUROBI_HOME"],"bin"), libgurobi)
