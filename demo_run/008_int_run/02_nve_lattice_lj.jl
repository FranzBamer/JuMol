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

## plot properties
size1 = 13.5
size2 = 25.0
hfig = 500
bfig = 680

folder_input = "demo_run/999_benchmark_samples/LJ_glass/10_10/"
folder_output = "demo_run/008_int_run/res/"


## create Jumol object
flush(stdout)
molstruc = Jumol.Structure()
molstruc.rc = 2.4
molstruc.rskin = 0.4
molstruc.pbx = 1
molstruc.pby = 1
#
Jumol.initialize_structure_objects(molstruc)

## create LJ glass
Generate_lat = Jumol.Gen_lattice(molstruc)
Jumol.create_LJ_lattice(Generate_lat, 6, 3)

## initial box dimensions
lx0 = molstruc.box.lx
ly0 = molstruc.box.ly

## move one atom (should go back during minimization)
molstruc.atom_list[13].pos[1] += 0.3

Jumol.initialize_structure(molstruc)
flush(stdout)
molcalc = Jumol.Calc(molstruc)
Jumol.initialize_potential(molcalc,33)

## assign masses to the particles
for i in 1:molstruc.noa
        molstruc.atom_list[i].mass = 1.0
end

## set units
molunits = Jumol.Units(0)
Jumol.set_units(molunits)

## plot atomic structure
Plotter = Jumol.Vis2d(molstruc)
Plotter.bfig = bfig
Plotter.hfig = hfig
fig = Jumol.plot_box(Plotter)
Jumol.plot_atomic_structure_binary(Plotter, size1, size2,
                      [-lx0*0.075, lx0*1.075, lx0*1.075,-lx0*0.075],
                      [-ly0*0.075, -ly0*0.075, ly0*1.075, ly0*1.075] )
display(fig)

## time integration
clctr = Jumol.open_file(molstruc.reader,join([folder_output,"LJ_lattice_output.lammpstrj"]),"w")
Jumol.write_box(molstruc.reader,molstruc,1)
molint = Jumol.Integ(molstruc,molcalc,molunits,0.001)
#Jumol.run_nve(molint,1000,"nve_anim.gif",true)
Jumol.run_nve(molint,1000)
#Jumol.write_box(molstruc.reader,molstruc,2)
Jumol.close_file(molstruc.reader)

## Done
println("...Done...")
