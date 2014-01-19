using BinDeps

@BinDeps.setup

libgurobi = library_dependency("libgurobi",aliases=["libgurobi51.so","gurobi51","libgurobi55.so","gurobi55","libgurobi56.so","gurobi56"])

if haskey(ENV, "GUROBI_HOME")
    @unix_only provides(Binaries, joinpath(ENV["GUROBI_HOME"],"lib"), libgurobi)
end
@windows_only provides(Binaries, joinpath(ENV["GUROBI_HOME"],"bin"), libgurobi)

@BinDeps.install [:libgurobi => :libgurobi]
