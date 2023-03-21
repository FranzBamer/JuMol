## run Tersoff potential test file
######################################
######################################
## info
# this file executes only two atoms to test
# the Tersoff potential function, two atoms are initially
# positioned with a distance of 0.5 A and then the distance
# is increased by moving the second atom
######################################
######################################

Juno.clearconsole()

using PrettyTables
using Plots

include("../../../src/Jumol.jl")
flush(stdout)

molstruc = Jumol.Structure()
molstruc.rc = 2.1          # 2.1
molstruc.rskin = 1.0
molstruc.pbx = 0
molstruc.pby = 0

## input -- output
output_folder = "demo_run/002_potential_file_runs/Tersoff_C/"

#
Jumol.initialize_structure_objects(molstruc)
Jumol.create_box_by_hand(molstruc,20.0,20.0,1.0,0.0,0.0,0.0)
Jumol.add_atom_by_hand(molstruc,1,1,5.0,10.0,0.0,0.0,0.0,0.0)
Jumol.add_atom_by_hand(molstruc,2,1,6.0,10.0,0.0,0.0,0.0,0.0)
Jumol.add_atom_by_hand(molstruc,3,1,8.0,10.0,0.0,0.0,0.0,0.0)
#
Jumol.initialize_structure(molstruc)
#
molcalc = Jumol.Calc(molstruc)
Jumol.initialize_potential(molcalc,21)

## plot the structure
Plotter = Jumol.Vis2d(molstruc)
Plotter.bfig = 500
Plotter.hfig = 500
fig_sample = Jumol.plot_box(Plotter)
Jumol.plot_atoms(Plotter,"blue",5)
savefig(fig_sample,join([output_folder,"sample.pdf"]))
display(fig_sample)


## potential and force C-C
delta = 0.005
pot_vec_C_C = zeros(0)
force_vec_C_C = zeros(0)
r_vec_C_C = zeros(0)

for i in 1:200
	Jumol.calc_all_pair_forces_Tersoff(molcalc)
	molstruc.atom_list[2].pos[1] += delta
	Jumol.update_distances(molstruc)
	Jumol.calc_all_pair_pot(molcalc)
	Jumol.calc_all_pair_forces_Tersoff(molcalc)
	push!(r_vec_C_C,abs(molstruc.atom_list[2].pos[1]-molstruc.atom_list[1].pos[1]))
	push!(pot_vec_C_C,molstruc.atom_list[2].pot)
	push!(force_vec_C_C,molstruc.atom_list[2].force[1])
end
#
pot_C_C = plot(size=(500,400),xlabel="r",ylabel="pot",legend=:topright)
plot!(r_vec_C_C[1:end],pot_vec_C_C[1:end],color=:blue,linewidth=3,label="pot: C_C")
savefig(pot_C_C,join([output_folder,"pot_C_C.pdf"]))
display(pot_C_C)

force_C_C = plot(size=(500,400),xlabel="r",ylabel="force",legend=:topright)
plot!(r_vec_C_C[1:end],force_vec_C_C[1:end],color=:red,linewidth=3,label="force: C_C")
savefig(force_C_C,join([output_folder,"force_C_C.pdf"]))
display(force_C_C)

print("...Done...")
