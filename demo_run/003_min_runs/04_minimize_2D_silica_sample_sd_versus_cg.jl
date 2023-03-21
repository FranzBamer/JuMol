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
molstruc.rskin = 0.5
molstruc.pbx = 1
molstruc.pby = 1
Jumol.initialize_structure_objects(molstruc)

## input -- output
input_folder = "demo_run/999_benchmark_samples/2D_silica_8x8/sig_1.00/"

## plot properties
size1 = 2.75
size2 = 3.5
hfig = 500
bfig = 680

## load atomic structure
Jumol.read_lammpstrj(molstruc,join([input_folder,"sample_0.lammpstrj"]))

## initialize the structure
Jumol.initialize_structure(molstruc)
flush(stdout)
molcalc = Jumol.Calc(molstruc)
Jumol.initialize_potential(molcalc,13)

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
start = time()
Jumol.calc_all_pair_forces(molcalc)
Minimizer = Jumol.Min(molstruc,molcalc)
Minimizer.alpha_min = 1.0e0
Jumol.run_cg(Minimizer,500,1.0e-5)
#Jumol.run_cg(Minimizer,10)
elapsed = time() - start
println("time: ",elapsed)

## plot atomic structure again
Plotter = Jumol.Vis2d(molstruc)
Plotter.bfig = bfig
Plotter.hfig = hfig
fig = Jumol.plot_box(Plotter)
Jumol.plot_atomic_structure_binary(Plotter, size1, size2,
                      [-molstruc.box.lx*0.075, molstruc.box.lx*1.075, molstruc.box.lx*1.075,-molstruc.box.lx*0.075],
                      [-molstruc.box.ly*0.075, -molstruc.box.ly*0.075, molstruc.box.ly*1.075, molstruc.box.ly*1.075] )
display(fig)







println("...Done...")
