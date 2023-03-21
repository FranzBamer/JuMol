######################################
######################################
######################################
## set the box and the atoms to affine deformations
# 2D silica benchmark sample

Juno.clearconsole()



using PrettyTables
using Plots


include("../../src/Jumol.jl")
flush(stdout)

molstruc = Jumol.Structure()
molstruc.rc = 10.0
molstruc.rskin = 0.5
molstruc.pbx = 1
molstruc.pby = 1

## plot properties
size1 = 2.75
size2 = 3.5
hfig = 500
bfig = 680

## input -- output
input_folder = "demo_run/999_benchmark_samples/2D_silica_8x8/sig_1.00/"

Jumol.initialize_structure_objects(molstruc)
filename = join([input_folder,"sample_0.lammpstrj"])
Jumol.read_lammpstrj(molstruc,filename)
Jumol.initialize_structure(molstruc)

lx0 = molstruc.box.lx
ly0 = molstruc.box.ly

## plot atomic structure deformed
Plotter = Jumol.Vis2d(molstruc)
Plotter.bfig = bfig
Plotter.hfig = hfig
fig = Jumol.plot_box(Plotter)
Jumol.plot_atomic_structure_binary(Plotter, size1, size2,
                      [-lx0*0.075, lx0*1.50, lx0*1.50,-lx0*0.075],
                      [-ly0*0.075, -ly0*0.075, ly0*1.075, ly0*1.075] )
display(fig)

Deformer = Jumol.Aff_deform(molstruc)
Jumol.set_affine_deform_shear(Deformer,20.0,0.0,0.0)

## plot atomic structure deformed
Plotter = Jumol.Vis2d(molstruc)
Plotter.bfig = bfig
Plotter.hfig = hfig
fig = Jumol.plot_box(Plotter)
Jumol.plot_atomic_structure_binary(Plotter, size1, size2,
                      [-lx0*0.075, lx0*1.50, lx0*1.50,-lx0*0.075],
                      [-ly0*0.075, -ly0*0.075, ly0*1.075, ly0*1.075] )
display(fig)

## Done
println("...Done...")
