## run a minimization of one configuration
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
folder_input = "demo_run/999_benchmark_samples/LJ_glass/10_10/"
output_file = "demo_run/006_AQS_run/03_tension/res/"

## plot properties (for LJ glass 10x10 to be chosen appropriately)
size1 = 13.5
size2 = 25.0
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

## load LJ glass
Jumol.read_lammpstrj(molstruc,join([folder_input,"sample_10_10_num_2.lammpstrj"]))

## initialize the structure
Jumol.initialize_structure(molstruc)
flush(stdout)
molcalc = Jumol.Calc(molstruc)
Jumol.initialize_potential(molcalc,33)

## initial box dimensions
lx0 = molstruc.box.lx
ly0 = molstruc.box.ly

## plot atomic structure
Plotter = Jumol.Vis2d(molstruc)
Plotter.bfig = bfig
Plotter.hfig = hfig
fig = Jumol.plot_box(Plotter)
Jumol.plot_atomic_structure_binary(Plotter, size1, size2,
                      [-lx0*0.075, lx0*1.40, lx0*1.40,-lx0*0.075],
                      [-ly0*0.075, -ly0*0.075, ly0*1.075, ly0*1.075] )
display(fig)

## initial minimizer
alpha_min = 1.0e0
#
time_start = time()
eps = 1.0e-7
max_num_steps = 1000000
Minimizer = Jumol.Min(molstruc,molcalc)
Minimizer.alpha_min = alpha_min
Minimizer.acc_factor = 1.001
Jumol.run_cg(Minimizer,max_num_steps,eps)
lx0 = molstruc.box.lx
ly0 = molstruc.box.ly
time_end = time()
println("time elapsed: ", time_end-time_start)

## deformation object
Deformer = Jumol.Aff_deform(molstruc)

## run AQS steps
time_start = time()
num_steps = 50
delta_x = 0.01
stress_response_mat = zeros(num_steps,2)
fname = join([output_file, "output_hist_tension.lammpstrj"])
Jumol.open_file(molstruc.reader, fname, "w")
stress_strain_file = open(join([output_file, "stress_strain_tension.dat"]), "w")
anim = @animate for i in 1:num_steps
    println("AQS step: ", i)
    ## deform
    Jumol.set_affine_deform_vol(Deformer,delta_x,0.0,0.0)
    # create fracture position
    molstruc.atom_list[55].pos[1] += 1.0e-5
    ## update distance
    Jumol.update_distances(molstruc)
    Jumol.update_neighbor_lists(molstruc)
    ## minimize
    Minimizer.alpha_min = alpha_min
    Minimizer.acc_factor = 1.001
    Jumol.run_cg(Minimizer,max_num_steps,eps)
    Jumol.put_atoms_back_to_box(molstruc)
    ## plot atomic structure
    Plotter = Jumol.Vis2d(molstruc)
    Plotter.bfig = bfig
    Plotter.hfig = hfig
    fig = Jumol.plot_box(Plotter)
    Jumol.plot_atomic_structure_binary(Plotter, size1, size2,
                          [-lx0*0.075, lx0*1.40, lx0*1.40,-lx0*0.075],
                          [-ly0*0.075, -ly0*0.075, ly0*1.075, ly0*1.075] )
    display(fig)
    ## write box output
    Jumol.write_box(molstruc.reader, molstruc, i)
    flush(molstruc.reader.output_file)
    ## calculate stresses
    Jumol.calc_all_pair_forces(molcalc,true)
    println("stress tensor:")
    pretty_table(molcalc.stress_tensor, noheader = true,
        crop = :horizontal, formatters = ft_round(8))
    stress_response_mat[i,1] = molstruc.box.lx - lx0
    stress_response_mat[i,2] = molcalc.stress_tensor[1,1]
    write(stress_strain_file, string(stress_response_mat[i,1]," ", string(stress_response_mat[i,2],"\n")))
    flush(stress_strain_file)
end
Jumol.close_file(molstruc.reader)
close(stress_strain_file)
elapsed_time = time() - time_start
println("calculation time: ", elapsed_time)
gif(anim, join([output_file,"video_deform_tension.gif"]), fps = 25)

## plot stress-strain result
fig_stress_strain = plot(stress_response_mat[:,1],stress_response_mat[:,2],color="blue",lw=1,fmt=:pdf)
display(fig_stress_strain)

##
println("...Done...")
