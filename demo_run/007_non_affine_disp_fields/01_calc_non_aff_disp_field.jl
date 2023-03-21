## run a minimization of one configuration
######################################
######################################
######################################

Juno.clearconsole()


using PrettyTables
using Plots
using DelimitedFiles

## load module Jumol
include("../../src/Jumol.jl")
flush(stdout)

## input -- output
folder_input = "demo_run/006_AQS_run/02_simple_shear/res/"
output_file = "demo_run/007_non_affine_disp_fields/res/"

##
init_step = 1
step1 = 10
step2 = 11
#
factor = 100.0
factor_lw = 60.0
limit_disp = 5.0

## plot properties (for LJ glass 10x10 to be chosen appropriately)
size1 = 13.5
size2 = 25.0
hfig = 500
bfig = 680

## create init Jumol object
flush(stdout)
ms_init = Jumol.Structure()
ms_init.rc = 2.4
ms_init.rskin = 0.4
ms_init.pbx = 1
ms_init.pby = 1
Jumol.initialize_structure_objects(ms_init)
## load LJ glass
Jumol.read_lammpstrj(ms_init,join([folder_input,"output_hist_LJ_glass.lammpstrj"]),init_step)
## initialize the structure
Jumol.initialize_structure(ms_init)
flush(stdout)

## initial box dimensions
lx0 = molstruc.box.lx
ly0 = molstruc.box.ly

## plot atomic structure
Plotter = Jumol.Vis2d(molstruc)
Plotter.bfig = bfig
Plotter.hfig = hfig
fig = Jumol.plot_box(Plotter)
Jumol.plot_atomic_structure_binary(Plotter, size1, size2,
                      [-lx0*0.075, lx0*1.3, lx0*1.3,-lx0*0.075],
                      [-ly0*0.075, -ly0*0.075, ly0*1.075, ly0*1.075] )
display(fig)


#### create step1 jumol object
flush(stdout)
ms_1 = Jumol.Structure()
ms_1.rc = 2.4
ms_1.rskin = 0.4
ms_1.pbx = 1
ms_1.pby = 1
Jumol.initialize_structure_objects(ms_1)
## load LJ glass
Jumol.read_lammpstrj(ms_1,join([folder_input,"output_hist_LJ_glass.lammpstrj"]),step1)
## initialize the structure
Jumol.initialize_structure(ms_1)
flush(stdout)
## plot atomic structure
Plotter = Jumol.Vis2d(molstruc)
Plotter.bfig = bfig
Plotter.hfig = hfig
fig = Jumol.plot_box(Plotter)
Jumol.plot_atomic_structure_binary(Plotter, size1, size2,
                      [-lx0*0.075, lx0*1.3, lx0*1.3,-lx0*0.075],
                      [-ly0*0.075, -ly0*0.075, ly0*1.075, ly0*1.075] )
display(fig)

## create step2 jumol object
flush(stdout)
ms_2 = Jumol.Structure()
ms_2.rc = 2.4
ms_2.rskin = 0.4
ms_2.pbx = 1
ms_2.pby = 1
Jumol.initialize_structure_objects(ms_2)
## load LJ glass
Jumol.read_lammpstrj(ms_2,join([folder_input,"output_hist_LJ_glass.lammpstrj"]),step2)
## initialize the structure
Jumol.initialize_structure(ms_2)
flush(stdout)
## plot atomic structure
Plotter = Jumol.Vis2d(molstruc)
Plotter.bfig = bfig
Plotter.hfig = hfig
fig = Jumol.plot_box(Plotter)
Jumol.plot_atomic_structure_binary(Plotter, size1, size2,
                      [-lx0*0.075, lx0*1.3, lx0*1.3,-lx0*0.075],
                      [-ly0*0.075, -ly0*0.075, ly0*1.075, ly0*1.075] )
display(fig)

## calculate non-affine displacement field
Df = Jumol.Disp_field(ms_init, ms_1, ms_2)
Df.limit_val = limit_disp
Jumol.get_coord_matrices(Df)
aff_fields = Jumol.calc_aff_disp_field(Df)
full_field = Jumol.calc_full_disp_field(Df)
non_aff_field = full_field - aff_fields[1]
#
Plotter = Jumol.Vis2d(ms_init)
Plotter.bfig = 500
Plotter.hfig = 500
fig = Jumol.plot_box(Plotter,"black")
Jumol.plot_displacement_field(Df, non_aff_field, aff_fields[2], factor, factor_lw)
savefig(join([output_file,"stz_",string(step1),"_",string(step2),".svg"]))
display(fig)
#


##
println("...Done...")
