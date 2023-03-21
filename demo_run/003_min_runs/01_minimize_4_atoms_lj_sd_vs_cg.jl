## run a minimization of one configuration
######################################
######################################
######################################

Juno.clearconsole()


#using Revise
using PrettyTables
using Plots

## load module Jumol
include("../../src/Jumol.jl")
flush(stdout)

## create Jumol object
flush(stdout)
molstruc = Jumol.Structure()
molstruc.rc = 2.5
molstruc.rskin = 0.4
molstruc.pbx = 1
molstruc.pby = 1
#
Jumol.initialize_structure_objects(molstruc)

## build atomic structure
Jumol.create_box_by_hand(molstruc,10.0,10.0,10.0,0.0,0.0,0.0)
#
Jumol.add_atom_by_hand(molstruc,1,1,1.0,1.0,0.0,0.0,0.0,0.0)
Jumol.add_atom_by_hand(molstruc,2,2,2.5,1.3,0.0,0.0,0.0,0.0)
Jumol.add_atom_by_hand(molstruc,3,1,2.0,3.0,0.0,0.0,0.0,0.0)
Jumol.add_atom_by_hand(molstruc,4,2,4.0,3.0,0.0,0.0,0.0,0.0)
#

Jumol.initialize_structure(molstruc)
flush(stdout)
molcalc = Jumol.Calc(molstruc)
Jumol.initialize_potential(molcalc,33)

## assign masses to the particles
for i in 1:molstruc.noa
        molstruc.atom_list[i].mass = 28.0855
end

## plot atomic structure
Plotter = Jumol.Vis2d(molstruc)
Plotter.bfig = 500
Plotter.hfig = 500
fig = Jumol.plot_box(Plotter)
Jumol.plot_atoms(Plotter,"blue",5)
display(fig)

## run minimization either sd or cg
Minimizer = Jumol.Min(molstruc,molcalc)

## units
molunits = Jumol.Units(0)
Jumol.set_units(molunits)
#
#Jumol.run_sd(Minimizer,8000)
Minimizer.alpha_min = 1.0e0
Jumol.run_cg(Minimizer,500)


## plot atomic structure again
Plotter = Jumol.Vis2d(molstruc)
Plotter.bfig = 500
Plotter.hfig = 500
fig = Jumol.plot_box(Plotter)
Jumol.plot_atoms(Plotter,"blue",5)
display(fig)

## print result
println()
for i in 1:molstruc.noa
        println(molstruc.atom_list[i])
end







println("...Done...")
