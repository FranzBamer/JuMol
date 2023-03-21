## run a simple shear AQS run
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

## input -- output
output_file = "demo_run/006_AQS_run/02_simple_shear/res/"

## plot properties
size1 = 30.0
size2 = 40.0
hfig = 500
bfig = 680

## create Jumol object
flush(stdout)
molstruc = Jumol.Structure()
molstruc.rc = 2.4
molstruc.rskin = 0.4
molstruc.pbx = 1
molstruc.pby = 1
Jumol.initialize_structure_objects(molstruc)

## create LJ glass
Generate_lat = Jumol.Gen_lattice(molstruc)
Jumol.create_LJ_lattice(Generate_lat, 6, 3)

## initialize the structure
Jumol.initialize_structure(molstruc)
flush(stdout)
molcalc = Jumol.Calc(molstruc)
Jumol.initialize_potential(molcalc,33)

## plot atomic structure
Plotter = Jumol.Vis2d(molstruc)
Plotter.bfig = bfig
Plotter.hfig = hfig
fig = Jumol.plot_box(Plotter)
Jumol.plot_atomic_structure_binary(Plotter, size1, size2,
                      [-molstruc.box.lx*0.075, molstruc.box.lx*1.40, molstruc.box.lx*1.40,-molstruc.box.lx*0.075],
                      [-molstruc.box.ly*0.075, -molstruc.box.ly*0.075, molstruc.box.ly*1.075, molstruc.box.ly*1.075] )
display(fig)

## initial minimizer
alpha_min = 1.0e-1
time_start = time()
eps = 1.0e-7
max_num_steps = 100000
Minimizer = Jumol.Min(molstruc,molcalc)
Minimizer.alpha_min = alpha_min
Minimizer.acc_factor = 1.001
Jumol.run_cg(Minimizer,max_num_steps,eps)
lxy0 = molstruc.box.lxy
time_end = time()
println("time elapsed: ", time_end-time_start)

## deformation object
Deformer = Jumol.Aff_deform(molstruc)

## run AQS steps
time_start = time()
num_steps = 80
delta_lxy = 0.01
stress_response_mat = zeros(num_steps,2)
Jumol.open_file(molstruc.reader, join([output_file,"output_hist_cryst.lammpstrj"]), "w")
stress_strain_file = open(join([output_file,"stress_strain_output_LJ_cryst.dat"]), "w")
anim = @animate for i in 1:num_steps
    println("AQS step: ", i)
    ## deform
    Jumol.set_affine_deform_shear(Deformer,delta_lxy,0.0,0.0)
    # create fracture position
    molstruc.atom_list[16].pos[1] += 1.0e-4
    ## update distance
    Jumol.update_distances(molstruc)
    Jumol.update_neighbor_lists(molstruc)
    ## minimize
    Minimizer.alpha_min = alpha_min
    Jumol.run_cg(Minimizer,max_num_steps,eps)
    #Jumol.put_atoms_back_to_box(molstruc)
    ## plot atomic structure
    Plotter = Jumol.Vis2d(molstruc)
    Plotter.bfig = bfig
    Plotter.hfig = hfig
    fig = Jumol.plot_box(Plotter)
    Jumol.plot_atomic_structure_binary(Plotter, size1, size2,
                          [-molstruc.box.lx*0.075, molstruc.box.lx*1.40, molstruc.box.lx*1.40,-molstruc.box.lx*0.075],
                          [-molstruc.box.ly*0.075, -molstruc.box.ly*0.075, molstruc.box.ly*1.075, molstruc.box.ly*1.075] )
    display(fig)
    ## write box output
    Jumol.write_box(molstruc.reader, molstruc, i)
    flush(molstruc.reader.output_file)
    ## calculate stresses
    Jumol.calc_all_pair_forces(molcalc,true)
    println("stress tensor:")
    pretty_table(molcalc.stress_tensor, noheader = true,
        crop = :horizontal, formatters = ft_round(8))
    stress_response_mat[i,1] = molstruc.box.lxy - lxy0
    stress_response_mat[i,2] = molcalc.stress_tensor[1,2]
    write(stress_strain_file, string(stress_response_mat[i,1]," ", string(stress_response_mat[i,2],"\n")))
    flush(stress_strain_file)
end
Jumol.close_file(molstruc.reader)
close(stress_strain_file)
elapsed_time = time() - time_start
println("calculation time: ", elapsed_time)

## save video
gif(anim, join([output_file,"benchmark_video_2D_LJ_cryst.gif"]), fps = 25)

## plot stress-strain result
fig_stress_strain = plot(stress_response_mat[:,1],stress_response_mat[:,2],color="blue",lw=1,fmt=:pdf)
display(fig_stress_strain)

##
println("...Done...")
