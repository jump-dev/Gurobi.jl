using Gurobi

dir = pwd()
cd(dir)

# A simple mip
#
#   maximize z
#
#   s.t.  2 p1 + 5 p2 + 4 p3 - z = 0
#         0.5 p1 + 2 p2 + 1 p3 <= 6
#
#         p1 is integer: 0 <= p1
#         p1 is integer: 0 <= p2
#         p1 is integer: 0 <= p3
#         z is binary

## BUILD MODEL
#-------------

simple_model_env = Gurobi.Env()

simple_model = Gurobi.Model(simple_model_env, "simple_mip", :maximize)

add_ivar!(simple_model, 0., 0, Inf)  # p1
add_ivar!(simple_model, 0., 0, Inf)  # p2
add_ivar!(simple_model, 0., 0, Inf)  # p3
add_cvar!(simple_model, 1., 0., Inf) # z
update_model!(simple_model)

# two constraints due to '='
add_constr!(simple_model, [2., 5., 4., -1.], '<', 0.)
add_constr!(simple_model, [2., 5., 4., -1.], '>', 0.)

add_constr!(simple_model, [0.5, 2., 1., 0.], '<', 6.)

setparam!(simple_model, "Heuristics", 0.0)
setparam!(simple_model, "Presolve", 0)

update_model!(simple_model)

optimize(simple_model)
println("-----------------------------")
println()

## WRITE SOLUTION FILE
#---------------------

write_model(simple_model, "simple_out.sol")
write_model(simple_model, "simple_out.mst")

## READ SOLUTION FILE
#--------------------

start_1 = Gurobi.get_dblattrarray(simple_model, "Start", 1, 4)
println("Current start vector:")
println(start_1)
println()

println("----------------------")
println("Reading solution file.")
println("----------------------")
println()
Gurobi.read(simple_model, "simple_out.sol")
update_model!(simple_model)

println("Updated start vector:")
start_2 = Gurobi.get_dblattrarray(simple_model, "Start", 1, 4)
println(start_2)

## VERIFY MIP START
#------------------

#The parameters Presolve and Heuristics must be set to 0
println()
println("-----------------------------")
optimize(simple_model)

## DELETE FILES
#--------------

rm("simple_out.sol")
rm("simple_out.mst")
