## run first system
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

include("../../../src/Jumol.jl")
flush(stdout)

## set current path
Jumol.set_the_current_path()

## distance vector
N = 1000
r_vec = zeros(N)
r_start = 0.6; r_end = 10.0
delta_r = (r_end-r_start)/float(N-1)
for i in 1:N
	r_vec[i] = r_start + (i-1)*delta_r
end


## structure of two atoms
molstruc = Jumol.Structure()
molstruc.rc = 10.0
molstruc.rskin = 1.0
molstruc.pbx = 0
molstruc.pby = 0
#
Jumol.initialize_structure_objects(molstruc)
Jumol.create_box_by_hand(molstruc,20.0,20.0,1.0,0.0,0.0,0.0)
Jumol.add_atom_by_hand(molstruc,1,1,5.0,10.0,0.0,0.0,0.0,0.0)
Jumol.add_atom_by_hand(molstruc,1,1,5.5,10.0,0.0,0.0,0.0,0.0)
#
Jumol.initialize_structure(molstruc)
#
molcalc = Jumol.Calc(molstruc)
Jumol.initialize_potential(molcalc,13)

## plot the structure
Plotter = Jumol.Vis2d(molstruc)
Plotter.bfig = 500
Plotter.hfig = 500
fig_sample = Jumol.plot_box(Plotter)
Jumol.plot_atoms(Plotter,"blue",5)
savefig(fig_sample,join(["sample.pdf"]))
display(fig_sample)


## potential and force Si-Si
delta = 0.01
pot_vec_Si_Si = zeros(0)
force_vec_Si_Si = zeros(0)
r_vec_Si_Si = zeros(0)
for i in 1:N
	Jumol.calc_all_pair_forces(molcalc)
	molstruc.atom_list[2].pos[1] = molstruc.atom_list[1].pos[1] + r_vec[i]
	Jumol.update_distances(molstruc)
	Jumol.calc_all_pair_pot(molcalc)
	Jumol.calc_all_pair_forces(molcalc)
	push!(r_vec_Si_Si,abs(molstruc.atom_list[2].pos[1]-molstruc.atom_list[1].pos[1]))
	push!(pot_vec_Si_Si,molstruc.atom_list[2].pot)
	push!(force_vec_Si_Si,molstruc.atom_list[2].force[1])
end
#
fig_Si_Si = plot(size=(500,400),xlabel="r",ylabel="pot / force")
plot!(r_vec_Si_Si[150:end],pot_vec_Si_Si[150:end],color=:blue,linewidth=3,label="pot Si-Si")
plot!(r_vec_Si_Si[175:end],force_vec_Si_Si[175:end],color=:red,linewidth=3,label="force: Si-Si")
savefig(fig_Si_Si,join("interaction_Si_Si.pdf"))
display(fig_Si_Si)


## potential and force Si-O
# resetting inital atom positions
molstruc.atom_list[2].pos[1] = 5.5
molstruc.atom_list[2].type = 2
pot_vec_Si_O = zeros(0)
force_vec_Si_O = zeros(0)
r_vec_Si_O = zeros(0)
for i in 1:N
	Jumol.calc_all_pair_forces(molcalc)
	molstruc.atom_list[2].pos[1] = molstruc.atom_list[1].pos[1] + r_vec[i]
	Jumol.update_distances(molstruc)
	Jumol.calc_all_pair_pot(molcalc)
	Jumol.calc_all_pair_forces(molcalc)
	push!(r_vec_Si_O,abs(molstruc.atom_list[2].pos[1]-molstruc.atom_list[1].pos[1]))
	push!(pot_vec_Si_O,molstruc.atom_list[2].pot)
	push!(force_vec_Si_O,molstruc.atom_list[2].force[1])
end
#
fig_Si_O = plot(size=(500,400),xlabel="r", ylabel="pot / force")
plot!(r_vec_Si_O[50:end], pot_vec_Si_O[50:end], color=:blue, linewidth=3, label="pot Si-O")
plot!(r_vec_Si_O[60:end], force_vec_Si_O[60:end], color=:red, linewidth=3, label="force: Si-O")
savefig(fig_Si_O,"interaction_Si_O.pdf")
display(fig_Si_O)


## potential and force O-O
# resetting inital atom positions
molstruc.atom_list[2].pos[1] = 5.5
molstruc.atom_list[1].type = 2
molstruc.atom_list[2].type = 2
delta = 0.01
pot_vec_O_O = zeros(0)
force_vec_O_O = zeros(0)
r_vec_O_O = zeros(0)
for i in 1:N
	Jumol.calc_all_pair_forces(molcalc)
	molstruc.atom_list[2].pos[1] = molstruc.atom_list[1].pos[1] + r_vec[i]
	Jumol.update_distances(molstruc)
	Jumol.calc_all_pair_pot(molcalc)
	Jumol.calc_all_pair_forces(molcalc)
	push!(r_vec_O_O,abs(molstruc.atom_list[2].pos[1]-molstruc.atom_list[1].pos[1]))
	push!(pot_vec_O_O,molstruc.atom_list[2].pot)
	push!(force_vec_O_O,molstruc.atom_list[2].force[1])
end
#
fig_O_O = plot(size=(500,400),xlabel="r",ylabel="pot / force")
plot!(r_vec_O_O[50:end], pot_vec_O_O[50:end], color=:blue, linewidth=3, label="pot O-O")
plot!(r_vec_O_O[60:end], force_vec_O_O[60:end], color=:red, linewidth=3, label="force: O-O")
savefig(fig_O_O, "interaction_O_O.pdf")
display(fig_O_O)



println("...Done...")
