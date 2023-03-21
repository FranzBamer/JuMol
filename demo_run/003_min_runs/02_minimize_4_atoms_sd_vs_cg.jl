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
molstruc.rc = 10.0
molstruc.rskin = 1.0
molstruc.pbx = 1
molstruc.pby = 1
Jumol.initialize_structure_objects(molstruc)

## build atomic structure
Jumol.create_box_by_hand(molstruc,12.0,12.0,10.0,0.0,0.0,0.0)
Jumol.add_atom_by_hand(molstruc,1,1,1.0,1.0,0.0,0.0,0.0,0.0)
Jumol.add_atom_by_hand(molstruc,2,2,4.0,5.0,0.0,0.0,0.0,0.0)
Jumol.add_atom_by_hand(molstruc,3,1,7.0,1.0,0.0,0.0,0.0,0.0)
Jumol.add_atom_by_hand(molstruc,4,1,9.0,8.0,0.0,0.0,0.0,0.0)

## initialize the structure
Jumol.initialize_structure(molstruc)
flush(stdout)
molcalc = Jumol.Calc(molstruc)
Jumol.initialize_potential(molcalc,13)

## plot atomic structure
Plotter = Jumol.Vis2d(molstruc)
Plotter.bfig = 500
Plotter.hfig = 500
fig = Jumol.plot_box(Plotter)
Jumol.plot_atoms(Plotter,"blue",5)
display(fig)

## run minimization (either sd or cg)
Minimizer = Jumol.Min(molstruc,molcalc)
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







println("...Done...")
