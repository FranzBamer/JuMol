## run a simple shear calculation of a 2D network glass
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
molstruc.linked_cells_bool = false
Jumol.initialize_structure_objects(molstruc)

## input -- output
folder_input = "demo_run/999_benchmark_samples/2D_silica_8x8/sig_1.00/"
folder_output = "demo_run/006_AQS_run/02_simple_shear/res/"

## load atomic structure
Jumol.read_lammpstrj(molstruc,join([folder_input,"sample_1.lammpstrj"]))

## initialize the structure
Jumol.initialize_structure(molstruc)
flush(stdout)
molcalc = Jumol.Calc(molstruc)
Jumol.initialize_potential(molcalc,13)

## plot properties
size1 = 2.75
size2 = 3.5
hfig = 500
bfig = 750

## plot atomic structure before minimization
Plotter = Jumol.Vis2d(molstruc)
Plotter.bfig = bfig
Plotter.hfig = hfig
fig = Jumol.plot_box(Plotter)
Jumol.plot_atomic_structure_binary(Plotter, size1, size2,
                      [-molstruc.box.lx*0.075, molstruc.box.lx*1.40, molstruc.box.lx*1.40,-molstruc.box.lx*0.075],
                      [-molstruc.box.ly*0.075, -molstruc.box.ly*0.075, molstruc.box.ly*1.075, molstruc.box.ly*1.075] )
display(fig)

## initial minimizer
time_start = time()
eps = 1.0e-5
max_num_steps = 10000
Minimizer = Jumol.Min(molstruc,molcalc)
Minimizer.alpha_min = 1.0e0
Minimizer.acc_factor = 1.001
Jumol.run_cg(Minimizer,max_num_steps,eps)
time_end = time()
println("time elapsed: ", time_end-time_start)

## plot atomic structure after minimization
Plotter = Jumol.Vis2d(molstruc)
Plotter.bfig = bfig
Plotter.hfig = hfig
fig = Jumol.plot_box(Plotter)
Jumol.plot_atomic_structure_binary(Plotter, size1, size2,
                      [-molstruc.box.lx*0.075, molstruc.box.lx*1.40, molstruc.box.lx*1.40,-molstruc.box.lx*0.075],
                      [-molstruc.box.ly*0.075, -molstruc.box.ly*0.075, molstruc.box.ly*1.075, molstruc.box.ly*1.075] )
display(fig)

## deformation object
Deformer = Jumol.Aff_deform(molstruc)

## run AQS steps
time_start = time()
num_steps = 10
delta_lxy = 0.1
#delta_x = 0.1
stress_response_mat = zeros(num_steps,2)
anim = @animate for i in 1:num_steps
    println("AQS step: ", i)
    ## deform
    Jumol.set_affine_deform_shear(Deformer,delta_lxy,0.0,0.0)
    #Jumol.set_affine_deform_vol(Deformer,delta_x,0.0,0.0)
    ## update distance
    Jumol.update_distances(molstruc)
    Jumol.update_neighbor_lists(molstruc)
    ## minimize
    Minimizer.alpha_min = 1.0e0
    Minimizer.acc_factor = 1.001
    Jumol.run_cg(Minimizer,max_num_steps,eps)
    ## plot atomic structure
    Plotter = Jumol.Vis2d(molstruc)
    Plotter.bfig = bfig
    Plotter.hfig = hfig
    fig = Jumol.plot_box(Plotter)
    Jumol.plot_atomic_structure_binary(Plotter, size1, size2,
                          [-molstruc.box.lx*0.075, molstruc.box.lx*1.40, molstruc.box.lx*1.40,-molstruc.box.lx*0.075],
                          [-molstruc.box.ly*0.075, -molstruc.box.ly*0.075, molstruc.box.ly*1.075, molstruc.box.ly*1.075] )
    display(fig)
    ## calculate stresses
    Jumol.calc_all_pair_forces(molcalc,true)
    println("stress tensor:")
    pretty_table(molcalc.stress_tensor, noheader = true,
        crop = :horizontal, formatters = ft_round(8))
    stress_response_mat[i,1] = molstruc.box.lxy
    stress_response_mat[i,2] = molcalc.stress_tensor[1,2]
end
elapsed_time = time() - time_start
println("calculation time: ", elapsed_time)
gif(anim, join([folder_output,"benchmark_video_2D_network_glass.gif"]), fps = 25)

## plot stress-strain result
fig_stress_strain = plot(stress_response_mat[:,1],stress_response_mat[:,2],color="blue",lw=1,fmt=:pdf)
#scatter!(stress_response_mat[:,1],stress_response_mat[:,2],color="blue")
display(fig_stress_strain)

## save stress_strain function
open(join([folder_output,"stress_strain_output_network_glass.dat"]),"w") do io
    writedlm(io,stress_response_mat)
end





println("...Done...")
