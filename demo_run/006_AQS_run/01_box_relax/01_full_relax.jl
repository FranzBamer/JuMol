## Example of a full relaxation of a box
######################################
######################################
######################################


Juno.clearconsole()



using PrettyTables
using Plots
using DelimitedFiles

## load module Jumol
include("../../../src/Jumol.jl")
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

## load atomic structure
Jumol.read_lammpstrj(molstruc,join([input_folder,"sample_2.lammpstrj"]))

## initialize the structure
Jumol.initialize_structure(molstruc)
flush(stdout)
molcalc = Jumol.Calc(molstruc)
Jumol.initialize_potential(molcalc,13)
lx0 = molstruc.box.lx
ly0 = molstruc.box.ly


## initial minimizer
epsilon = 1.0e-5
max_num_steps = 10000
Minimizer = Jumol.Min(molstruc,molcalc)
Minimizer.alpha_min = 1.0e0
Jumol.run_cg(Minimizer,max_num_steps,epsilon)

## plot initial stress tensor
Jumol.calc_all_pair_forces(molcalc,true)
println("initial stress_tensor:")
pretty_table(molcalc.stress_tensor, tf=tf_borderless,
             noheader = true, crop = :horizontal, formatters = ft_round(8))



## plot atomic structure initial
Plotter = Jumol.Vis2d(molstruc)
Plotter.bfig = bfig
Plotter.hfig = hfig
fig = Jumol.plot_box(Plotter)
Jumol.plot_atomic_structure_binary(Plotter, size1, size2,
                [-lx0*0.075, lx0*1.075, lx0*1.075,-lx0*0.075],
                [-ly0*0.075, -ly0*0.075, ly0*1.075, ly0*1.075] )
display(fig)

## deformation object
Deformer = Jumol.Aff_deform(molstruc)

## volumetric relaxation
## run AQS steps
num_steps = 1
delta_lxy = 0.1
stress_response_mat = zeros(num_steps,2)
anim = @animate for i in 1:num_steps
    println("###############")
    println("relaxation step: ", i)
    Borelax = Jumol.Box_relax(molcalc,Deformer,100,epsilon)
    Borelax.stress_tol = 1.0e-8
    Jumol.relax_biaxial(Borelax,1.0e-2,1.0e-2)
    ## update distance
    Jumol.update_distances(molstruc)
end

## plot atomic structure after relaxation
Plotter = Jumol.Vis2d(molstruc)
Plotter.bfig = bfig
Plotter.hfig = hfig
fig = Jumol.plot_box(Plotter)
Jumol.plot_atomic_structure_binary(Plotter, size1, size2,
                [-lx0*0.075, lx0*1.075, lx0*1.075,-lx0*0.075],
                [-ly0*0.075, -ly0*0.075, ly0*1.075, ly0*1.075] )
display(fig)
