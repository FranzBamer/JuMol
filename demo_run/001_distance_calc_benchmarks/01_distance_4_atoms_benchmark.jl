## run benchmark distances
######################################
######################################
######################################
## this file is to calculate the distance matrices and the
# verlet lists
# tested functions: update_distances, update_verlet_lists
# a hand calculation is added as pdf
####

Juno.clearconsole()

using PrettyTables
using Plots

include("../../src/Jumol.jl")
flush(stdout)

molstruc = Jumol.Structure()

molstruc.rc = 2.0
molstruc.rskin = 0.4
molstruc.pbx = 1
molstruc.pby = 1

Jumol.initialize_structure_objects(molstruc)
Jumol.create_box_by_hand(molstruc,5.0,5.0,10.0,0.0,0.0,0.0)
Jumol.add_atom_by_hand(molstruc,1,1,1.0,1.0,0.0,0.0,0.0,0.0)
Jumol.add_atom_by_hand(molstruc,2,1,3.0,2.0,0.0,0.0,0.0,0.0)
Jumol.add_atom_by_hand(molstruc,3,1,2.0,4.0,0.0,0.0,0.0,0.0)
Jumol.add_atom_by_hand(molstruc,4,1,4.0,4.0,0.0,0.0,0.0,0.0)
#
Jumol.initialize_structure(molstruc)

println("distances for atom 1:")
pretty_table(molstruc.atom_list[1].distances,
			 noheader = true, crop = :horizontal, formatters = ft_round(3))
println("distances for atom 2:")
pretty_table(molstruc.atom_list[2].distances,
			 noheader = true, crop = :horizontal, formatters = ft_round(3))
println("distances for atom 3:")
pretty_table(molstruc.atom_list[3].distances,
			 noheader = true, crop = :horizontal, formatters = ft_round(3))
println("neighbor lists for atom 1")
println(molstruc.atom_list[1].neighbor_indices)
println("neighbor lists for atom 2")
println(molstruc.atom_list[2].neighbor_indices)
println("neighbor lists for atom 3")
println(molstruc.atom_list[3].neighbor_indices)
println("neighbor lists for atom 4")
println(molstruc.atom_list[4].neighbor_indices)


Plotter = Jumol.Vis2d(molstruc)
Plotter.bfig = 500
Plotter.hfig = 500
fig = Jumol.plot_box(Plotter)
Jumol.plot_atoms(Plotter,"blue",5)
display(fig)

## Done
println("...Done...")
