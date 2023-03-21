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
molstruc.rc = 2.4
molstruc.rskin = 0.4
molstruc.pbx = 1
molstruc.pby = 1
#
Jumol.initialize_structure_objects(molstruc)

## create lj lattice
Gen_Lat = Jumol.Gen_lattice(molstruc)
Jumol.create_LJ_lattice(Gen_Lat,5,3)

## plot properties
size1 = 30.0
size2 = 40.0
hfig = 500
bfig = 680

## move one atom (should go back during minimization)
molstruc.atom_list[13].pos[1] += 0.1

Jumol.initialize_structure(molstruc)
flush(stdout)
molcalc = Jumol.Calc(molstruc)
Jumol.initialize_potential(molcalc,33)

## assign masses to the particles
for i in 1:molstruc.noa
        molstruc.atom_list[i].mass = 28.0855
end

## plot atomic structure initial
Plotter = Jumol.Vis2d(molstruc)
Plotter.bfig = bfig
Plotter.hfig = hfig
fig = Jumol.plot_box(Plotter)
Jumol.plot_atomic_structure_binary(Plotter, size1, size2,
                      [-molstruc.box.lx*0.075, molstruc.box.lx*1.075, molstruc.box.lx*1.075,-molstruc.box.lx*0.075],
                      [-molstruc.box.ly*0.075, -molstruc.box.ly*0.075, molstruc.box.ly*1.075, molstruc.box.ly*1.075] )
display(fig)

## run minimization
Minimizer = Jumol.Min(molstruc,molcalc)

## units
molunits = Jumol.Units(0)
Jumol.set_units(molunits)

#Jumol.run_sd(Minimizer,8000)
Minimizer.alpha_min = 1.0e0
Jumol.run_cg(Minimizer,500)


## plot atomic structure initial
Plotter = Jumol.Vis2d(molstruc)
Plotter.bfig = bfig
Plotter.hfig = hfig
fig = Jumol.plot_box(Plotter)
Jumol.plot_atomic_structure_binary(Plotter, size1, size2,
                      [-molstruc.box.lx*0.075, molstruc.box.lx*1.075, molstruc.box.lx*1.075,-molstruc.box.lx*0.075],
                      [-molstruc.box.ly*0.075, -molstruc.box.ly*0.075, molstruc.box.ly*1.075, molstruc.box.ly*1.075] )
display(fig)



println("...Done...")
