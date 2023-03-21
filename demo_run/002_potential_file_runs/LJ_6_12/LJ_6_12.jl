## run 6-12 LJ potential for different distances
######################################
######################################
## info
# this file executes only two atoms to test
# the yukawa potential function, two atoms are initially
# positioned with a distance of 0.5 A and then the distance
# is increased by moving the second atom
######################################
######################################

Juno.clearconsole()

using PrettyTables
using Plots
using DelimitedFiles


include("../../../src/Jumol.jl")
flush(stdout)

## write a file
function write_output_file(filename, x, y)
	output_file = open(filename,"w")
	for i in 1:size(x,1)
		write(output_file, string(x[i]," ", string(y[i],"\n")))
	end
	close(output_file)
end

## create Jumol object
molstruc = Jumol.Structure()
molstruc.rc = 5.0
molstruc.rskin = 0.4
molstruc.pbx = 0
molstruc.pby = 0

## input -- output
output_folder = "demo_run/002_potential_file_runs/LJ_6_12/"
#
Jumol.initialize_structure_objects(molstruc)
Jumol.create_box_by_hand(molstruc,20.0,20.0,1.0,0.0,0.0,0.0)
Jumol.add_atom_by_hand(molstruc,1,1,5.0,10.0,0.0,0.0,0.0,0.0)
Jumol.add_atom_by_hand(molstruc,1,1,5.5,10.0,0.0,0.0,0.0,0.0)
#
Jumol.initialize_structure(molstruc)
#
molcalc = Jumol.Calc(molstruc)
Jumol.initialize_potential(molcalc,33)

#
num_plot = 250
#

## plot the structure
Plotter = Jumol.Vis2d(molstruc)
Plotter.bfig = 500
Plotter.hfig = 500
fig_sample = Jumol.plot_box(Plotter)
Jumol.plot_atoms(Plotter,"blue",5)
savefig(fig_sample,join([output_folder,"sample.pdf"]))
display(fig_sample)

## units
molunits = Jumol.Units(1)
Jumol.set_units(molunits)

## potential and force particle_1 -- particle_2
delta = 0.01
pot_vec_1_1 = zeros(0)
force_vec_1_1 = zeros(0)
r_vec_1_1 = zeros(0)
for i in 1:num_plot
	Jumol.calc_all_pair_forces(molcalc)
	molstruc.atom_list[2].pos[1] += delta
	Jumol.update_distances(molstruc)
	Jumol.calc_all_pair_pot(molcalc)
	Jumol.calc_all_pair_forces(molcalc)
	push!(r_vec_1_1,abs(molstruc.atom_list[2].pos[1]-molstruc.atom_list[1].pos[1]))
	push!(pot_vec_1_1,molstruc.atom_list[2].pot)
	push!(force_vec_1_1,molstruc.atom_list[2].force[1])
end
#
num_start = 7
num_start_force = 17
fig_1_1 = plot(size=(500,400),xlabel="r",ylabel="pot / force")
plot!(r_vec_1_1[num_start:end],pot_vec_1_1[num_start:end],color=:blue,linewidth=3,label="pot 1-1")
plot!(r_vec_1_1[num_start_force:end],force_vec_1_1[num_start_force:end],color=:red,linewidth=3,label="force: 1-1")
savefig(fig_1_1,join([output_folder,"interaction_1_1.pdf"]))
display(fig_1_1)
# save file
fname = join([output_folder,"interaction_1_1.dat"])
fname_force = join([output_folder,"interaction_der_1_1.dat"])
write_output_file(fname, r_vec_1_1[num_start:end], pot_vec_1_1[num_start:end])
write_output_file(fname_force, r_vec_1_1[num_start:end], force_vec_1_1[num_start:end])



## potential and force Si-O
molstruc.atom_list[2].pos[1] = 5.5
molstruc.atom_list[2].type = 2
delta = 0.01
pot_vec_1_2 = zeros(0)
force_vec_1_2 = zeros(0)
r_vec_1_2 = zeros(0)
for i in 1:num_plot
	Jumol.calc_all_pair_forces(molcalc)
	molstruc.atom_list[2].pos[1] += delta
	Jumol.update_distances(molstruc)
	Jumol.calc_all_pair_pot(molcalc)
	Jumol.calc_all_pair_forces(molcalc)
	push!(r_vec_1_2,abs(molstruc.atom_list[2].pos[1]-molstruc.atom_list[1].pos[1]))
	push!(pot_vec_1_2,molstruc.atom_list[2].pot)
	push!(force_vec_1_2,molstruc.atom_list[2].force[1])
end
#
num_start = 42
num_start_force = 57
fig_1_2 = plot(size=(500,400),xlabel="r",ylabel="pot / force")
plot!(r_vec_1_2[num_start:end],pot_vec_1_2[num_start:end],color=:blue,linewidth=3,label="pot 1-2")
plot!(r_vec_1_2[num_start_force:end],force_vec_1_2[num_start_force:end],color=:red,linewidth=3,label="force: 1-2")
savefig(fig_1_2,join([output_folder,"interaction_1_2.pdf"]))
display(fig_1_2)
# save file
fname = join([output_folder,"interaction_1_2.dat"])
fname_force = join([output_folder,"interaction_der_1_2.dat"])
write_output_file(fname, r_vec_1_2[num_start:end], pot_vec_1_2[num_start:end])
write_output_file(fname_force, r_vec_1_2[num_start:end], force_vec_1_2[num_start:end])

## potential and force O-O
molstruc.atom_list[2].pos[1] = 5.5
molstruc.atom_list[1].type = 2
molstruc.atom_list[2].type = 2
delta = 0.01
pot_vec_2_2 = zeros(0)
force_vec_2_2 = zeros(0)
r_vec_2_2 = zeros(0)
for i in 1:num_plot
	Jumol.calc_all_pair_forces(molcalc)
	molstruc.atom_list[2].pos[1] += delta
	Jumol.update_distances(molstruc)
	Jumol.calc_all_pair_pot(molcalc)
	Jumol.calc_all_pair_forces(molcalc)
	push!(r_vec_2_2,abs(molstruc.atom_list[2].pos[1]-molstruc.atom_list[1].pos[1]))
	push!(pot_vec_2_2,molstruc.atom_list[2].pot)
	push!(force_vec_2_2,molstruc.atom_list[2].force[1])
end
#
num_start = 53
num_start_force = 72
fig_2_2 = plot(size=(500,400),xlabel="r",ylabel="pot / force")
plot!(r_vec_2_2[num_start:end],pot_vec_2_2[num_start:end],color=:blue,linewidth=3,label="pot 2-2")
plot!(r_vec_2_2[num_start_force:end],force_vec_2_2[num_start_force:end],color=:red,linewidth=3,label="force: 2-2")
savefig(fig_2_2,join([output_folder,"interaction_2_2.pdf"]))
display(fig_2_2)
# save file
fname = join([output_folder,"interaction_2_2.dat"])
fname_force = join([output_folder,"interaction_der_2_2.dat"])
write_output_file(fname, r_vec_2_2[num_start:end], pot_vec_2_2[num_start:end])
write_output_file(fname_force, r_vec_2_2[num_start:end], force_vec_2_2[num_start:end])



println("...Done...")
